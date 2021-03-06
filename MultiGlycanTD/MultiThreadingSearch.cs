using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.engine.score;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using SpectrumData;
using SpectrumProcess.algorithm;
using SpectrumProcess.deisotoping;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MultiGlycanTD
{
    using GlycanFragments = Dictionary<FragmentType, List<string>>;
    public class MultiThreadingSearch
    {
        Counter readingCounter;
        Counter searchCounter;
        GlycanJson glycanJson;
        CompdJson compdJson;

        ConcurrentQueue<SearchTask> tasks;
        ConcurrentQueue<SearchTask> decoyTasks;
        ConcurrentDictionary<int, ISpectrum> tandemSpectra;
        ConcurrentDictionary<int, ISpectrum> decoyTandemSpectra;

        string msPath;
        List<SearchResult> targets = new List<SearchResult>();
        List<SearchResult> decoys = new List<SearchResult>();
        private readonly object resultLock = new object();
        private readonly double searchRange = 1;
        int taskSize = 0;
        // It is not likely that a target spectrum has very high charge
        // for glycan or very few peaks for fragments.
        int minPeaks = 30;   // sequencable spectrum with min num peaks
        int minCharge = 2;
        // cluster
        int K = 4;
        int maxIter = 1000;
        double difference = 0.01;
        // noise
        ConcurrentBag<double> intensityBag;
        double noise = 0;

        public MultiThreadingSearch(string msPath,
            Counter readingCounter, Counter searchCounter,
            GlycanJson glycanJson)
        {
            this.msPath = msPath;
            this.readingCounter = readingCounter;
            this.searchCounter = searchCounter;
            this.glycanJson = glycanJson;
            compdJson = glycanJson.Compound;

            // read spectrum
            tasks = new ConcurrentQueue<SearchTask>();
            decoyTasks = new ConcurrentQueue<SearchTask>();
            tandemSpectra = new ConcurrentDictionary<int, ISpectrum>();
            decoyTandemSpectra = new ConcurrentDictionary<int, ISpectrum>();
            GenerateTasks();
            GenerateDecoyTasks();
            taskSize = tasks.Count + decoyTasks.Count;
        }

        public List<SearchResult> Target()
        { return targets; }

        public List<SearchResult> Decoy()
        { return decoys; }

        void UpdateTask(List<SearchResult> t, List<SearchResult> d)
        {
            lock (resultLock)
            {
                targets.AddRange(t);
                decoys.AddRange(d);
            }
        }

        public void Run()
        {
            List<Task> searches = new List<Task>();
            for (int i = 0; i < SearchingParameters.Access.ThreadNums; i++)
            {
                Task LastTask = new Task(() => LocalSearch());
                LastTask.Start();
                searches.Add(LastTask);
            }
            Task.WaitAll(searches.ToArray());

            IGlycanScorer scorer;
            scorer = new GlycanScorerCluster(
                SearchingParameters.Access.ThreadNums,
                SearchingParameters.Access.Similarity,
                SearchingParameters.Access.BinWidth,
                K, maxIter, difference);

            scorer.Init(tandemSpectra, targets);
            scorer.Run();
            targets = scorer.Result();

            scorer.Init(decoyTandemSpectra, decoys);
            scorer.Run();
            decoys = scorer.Result();
        }

        public ConcurrentDictionary<int, ISpectrum> MSMSSpectra()
        {
            return tandemSpectra;
        }

        public ConcurrentDictionary<int, ISpectrum> DecoySpectra()
        {
            return decoyTandemSpectra;
        }

        void GenerateTasks()
        {
            MultiThreadingSearchHelper.GenerateSearchTasks(msPath, tasks, 
                tandemSpectra, readingCounter, minPeaks, SearchingParameters.Access.MaxCharge, minCharge, searchRange);
        }

        void GenerateDecoyTasks()
        {
            MultiThreadingSearchHelper.GenerateSearchTasks(SearchingParameters.Access.DecoyFile,
                    decoyTasks, decoyTandemSpectra, readingCounter, minPeaks, SearchingParameters.Access.MaxCharge, minCharge, searchRange);
        }

        void TaskLocalSearch(ref List<SearchResult> results,
            SearchTask task, GlycanPrecursorMatch precursorMatch,
            GlycanEnvelopeMatch envelopeMatch, IGlycanSearch glycanSearch, 
            SearchMetaData searchInfo)
        {
            foreach (double ion in SearchingParameters.Access.Ions)
            {
                // precursor match
                List<string> candidates = precursorMatch.Match(task.PrecursorMZ, task.Charge, ion);
                if (candidates.Count > 0)
                {
                    // isotopic envelope
                    if (task.MSPeaks != null)
                    {
                        candidates = envelopeMatch.Match(candidates, task.MSPeaks,
                            task.PrecursorMZ, task.Charge, ion);
                    }
                    // spectrum search
                    List<SearchResult> searched = glycanSearch.Search(
                        candidates, task.Spectrum.GetPeaks(), task.Charge, ion);

                    if (searched.Count > 0)
                    {
                        searchInfo.Commit(searched, task.PrecursorMZ, task.Charge,
                            task.Spectrum.GetScanNum(), task.Spectrum.GetRetention());
                        results.AddRange(searched);
                    }
                }
            }
        }

        void LocalSearch()
        {
            List<SearchResult> tempResults = new List<SearchResult>();
            List<SearchResult> tempDecoyResults = new List<SearchResult>();

            ISearch<string> searcher = new BucketSearch<string>(
                SearchingParameters.Access.MS1ToleranceBy,
                SearchingParameters.Access.MS1Tolerance);
            GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, compdJson);

            EnvelopeProcessor envelopeProcessor = new EnvelopeProcessor(
                SearchingParameters.Access.MS1ToleranceBy,
                SearchingParameters.Access.MS1Tolerance,
                searchRange);
            GlycanEnvelopeMatch envelopeMatch = new GlycanEnvelopeMatch(
                envelopeProcessor, compdJson,
                SearchingParameters.Access.MS1ToleranceBy,
                SearchingParameters.Access.MS1Tolerance
                );

            ISearch<GlycanFragments> searcher2 = new BucketSearch<GlycanFragments>(
                SearchingParameters.Access.MS2ToleranceBy,
                SearchingParameters.Access.MSMSTolerance);

            Averagine averagine = new Averagine(AveragineType.PermethylatedGlycan);
            if (glycanJson.Derivation == DerivationType.Native)
            {
                averagine = new Averagine(AveragineType.Glycan);
            }
            AveragineDeisotoping deisotoping = new AveragineDeisotoping(averagine,
                ToleranceBy.Dalton, 0.1);
            IGlycanSearch glycanSearch
                = new GlycanSearchDeisotoping(searcher2, glycanJson, deisotoping);
            //IGlycanSearch glycanSearch
            //    = new GlycanSearch(searcher2, glycanJson);
            SearchMetaData searchInfo = new SearchMetaData();

            // targets
            while (tasks.TryDequeue(out SearchTask task))
            {
                TaskLocalSearch(ref tempResults, task,
                    precursorMatch, envelopeMatch, glycanSearch, searchInfo);

                searchCounter.Add(taskSize);
            }

            // decoys
            while (decoyTasks.TryDequeue(out SearchTask task))
            {
                TaskLocalSearch(ref tempDecoyResults, task,
                    precursorMatch, envelopeMatch, glycanSearch, searchInfo);

                searchCounter.Add(taskSize);
            }

            UpdateTask(tempResults, tempDecoyResults);
        }

    }
}
