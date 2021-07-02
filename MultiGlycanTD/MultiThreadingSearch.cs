using SpectrumProcess.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.engine.score;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumData.Spectrum;
using SpectrumProcess;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SpectrumProcess.deisotoping;

namespace MultiGlycanTD
{
    using GlycanFragments = Dictionary<FragmentTypes, List<string>>;
    public class MultiThreadingSearch
    {
        Counter readingCounter;
        Counter searchCounter;
        GlycanJson glycanJson;
        CompdJson compdJson;

        ConcurrentQueue<SearchTask> tasks;
        ConcurrentQueue<SearchTask> decoyTasks;
        ConcurrentDictionary<int, ISpectrum> tandemSpectra;

        string msPath;
        List<SearchResult> targets = new List<SearchResult>();
        List<SearchResult> decoys = new List<SearchResult>();
        private readonly object resultLock = new object();
        private readonly double searchRange = 1;
        int taskSize = 0;
        // It is not likely that a target spectrum has very high charge
        // for glycan or very few peaks for fragments.
        int minPeaks = 30;   // sequencable spectrum with min num peaks
        int minCharage = 2;
        int maxCharge = 4;   // max charge of spectrum to consider
        // Generate decoy by delta_M > d
        // for glycan delta_M < max_d, since glycan fragmetns
        // differ by monosaccradie can be very similar.
        ToleranceBy DistanceType = ToleranceBy.Dalton;   // decoy distance dalton or ppm
        int minDistance = 1;   // d, delta M > d
        int maxDistance = 50;  // max distance to consider
        int randomSeed = 2;   // deterministic results.

        public int Seed { get; set; } = 2;

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

            Dictionary<int, ISpectrum> spectra = new Dictionary<int, ISpectrum>(tandemSpectra);

            GlycanScorer scorer = new GlycanScorer(spectra, targets,
                SearchingParameters.Access.MS2ToleranceBy,
                SearchingParameters.Access.MSMSTolerance);
            scorer.Run();
            targets = scorer.Result();

            scorer = new GlycanScorer(spectra, decoys,
                SearchingParameters.Access.MS2ToleranceBy,
                SearchingParameters.Access.MSMSTolerance);
            scorer.Run();
            decoys = scorer.Result();
        }

        void GenerateTasks()
        {
            if (Path.GetExtension(msPath) == ".mgf")
            {
                MGFSpectrumReader mgfReader = new MGFSpectrumReader();
                mgfReader.Init(msPath);

                Dictionary<int, MS2Spectrum> spectraData = mgfReader.GetSpectrum();
                foreach (int scan in spectraData.Keys)
                {
                    MS2Spectrum spectrum = spectraData[scan];
                    if (spectrum.GetPeaks().Count <= minPeaks)
                        continue;

                    readingCounter.Add(spectraData.Count);

                    tandemSpectra[scan] = spectrum;
                    SearchTask searchTask = new SearchTask(spectrum,
                        spectrum.PrecursorMZ(), spectrum.PrecursorCharge());
                    tasks.Enqueue(searchTask);
                }
            }

            // read spectrum
            ISpectrumReader reader = new ThermoRawSpectrumReader();
            IProcess process = new WeightedAveraging(new LocalNeighborPicking());

            reader.Init(msPath);

            int start = reader.GetFirstScan();
            int end = reader.GetLastScan();

            Dictionary<int, List<int>> scanGroup = new Dictionary<int, List<int>>();
            int current = -1;
            for (int i = start; i < end; i++)
            {
                if (reader.GetMSnOrder(i) == 1)
                {
                    current = i;
                    scanGroup[i] = new List<int>();
                }
                else if (reader.GetMSnOrder(i) == 2
                    && reader.GetActivation(i) == TypeOfMSActivation.CID)
                {
                    scanGroup[current].Add(i);
                }
            }

            Parallel.ForEach(scanGroup,
                new ParallelOptions { MaxDegreeOfParallelism = SearchingParameters.Access.ThreadNums },
                (scanPair) =>
                {
                    if (scanPair.Value.Count > 0)
                    {
                        int scan = scanPair.Key;
                        ISpectrum ms1 = reader.GetSpectrum(scan);
                        if (ms1.GetPeaks().Count > 0)
                        {
                            foreach (int i in scanPair.Value)
                            {
                                // precurosr mz
                                double mz = reader.GetPrecursorMass(i, reader.GetMSnOrder(i));

                                // read ms1 peaks arouond precursor mz
                                List<IPeak> ms1Peaks = 
                                    MultiThreadingSearchHelper.FilterPeaks(ms1.GetPeaks(), mz, searchRange);
                                if (ms1Peaks.Count() == 0)
                                    continue;

                                // charage
                                ICharger charger = new Patterson();
                                int charge = charger.Charge(ms1Peaks, mz - searchRange, mz + searchRange);
                                if (charge > maxCharge || charge < minCharage)
                                    continue;

                                // search
                                ISpectrum ms2 = reader.GetSpectrum(i);
                                if (ms2.GetPeaks().Count <= minPeaks)
                                    continue;
                                ms2 = process.Process(ms2);
                                tandemSpectra[i] = ms2;
                                SearchTask searchTask = new SearchTask(ms2, mz, charge);
                                tasks.Enqueue(searchTask);
                                
                            }
                        }
                            
                    }
                    readingCounter.Add(scanGroup.Count);
                });
        }

        void GenerateDecoyTasks()
        {
            // find max precursor charges
            int maxCharge = 0;
            foreach (SearchTask task in tasks)
            {
                maxCharge = Math.Max(maxCharge, task.Charge);
            }

            // split spectrum on charges
            for (int charge = 0; charge <= maxCharge; charge++)
            {
                // init searcher to find all spectrum within delta < d
                ISearch<SearchTask> searcher = new BucketSearch<SearchTask>(DistanceType, maxDistance);
                List<Point<SearchTask>> points = new List<Point<SearchTask>>();
                foreach (SearchTask task in tasks)
                {
                    if (task.Charge == charge)
                    {
                        Point<SearchTask> point = new Point<SearchTask>(task.PrecursorMZ, task);
                        points.Add(point);
                    }
                }
                searcher.Init(points);

                // init randomness
                Random random = new Random(randomSeed);

                HashSet<int> swappedScans = new HashSet<int>();
                // precursor swap, swap any spectrums within d, set precursor mz
                foreach (SearchTask task in tasks)
                {
                    if (task.Charge != charge)
                        continue;

                    // avoid duplicate
                    if (swappedScans.Contains(task.Spectrum.GetScanNum()))
                        continue;

                    List<SearchTask> candidates = searcher.SearchContent(task.PrecursorMZ);
                    foreach (SearchTask selectTask in candidates)
                    {
                        // avoid duplicate
                        if (swappedScans.Contains(selectTask.Spectrum.GetScanNum()))
                            continue;

                        // distance bound by minDistance
                        if (DistanceType == ToleranceBy.PPM)
                        {
                            if (Math.Abs(selectTask.PrecursorMZ - task.PrecursorMZ)
                                / task.PrecursorMZ * 1000000.0 <= minDistance)
                                continue;
                        }
                        else
                        {
                            if (Math.Abs(selectTask.PrecursorMZ - task.PrecursorMZ) < minDistance)
                                continue;
                        }

                        // random picked
                        int r = random.Next(0, 2);
                        if (r % 2 == 0) continue;

                        // swap
                        decoyTasks.Enqueue(new SearchTask(task.Spectrum, selectTask.PrecursorMZ, charge));
                        decoyTasks.Enqueue(new SearchTask(selectTask.Spectrum, task.PrecursorMZ, charge));
                        swappedScans.Add(task.Spectrum.GetScanNum());
                        swappedScans.Add(selectTask.Spectrum.GetScanNum());
                        break;
                    }
                }
            };
        }

        void TaskLocalSearch(ref List<SearchResult> results,
            SearchTask task, GlycanPrecursorMatch precursorMatch,
            IGlycanSearch glycanSearch, SearchMetaData searchInfo)
        {
            foreach (double ion in SearchingParameters.Access.Ions)
            {
                //precursor match
                List<string> candidates = precursorMatch.Match(task.PrecursorMZ, task.Charge, ion);
                if (candidates.Count > 0)
                {
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
          

            ISearch<GlycanFragments> searcher2 = new BucketSearch<GlycanFragments>(
                SearchingParameters.Access.MS2ToleranceBy, 
                SearchingParameters.Access.MSMSTolerance);

            Averagine averagine = new Averagine(AveragineType.PermethylatedGlycan);
            if (glycanJson.Derivation == DerivationType.Native)
            {
                averagine = new Averagine(AveragineType.Glycan);
            }
            AveragineDeisotoping deisotoping = new AveragineDeisotoping(averagine,
                maxCharge, SearchingParameters.Access.MS2ToleranceBy,
                SearchingParameters.Access.MSMSTolerance);
            IGlycanSearch glycanSearch 
                = new GlycanSearchDeisotoping(searcher2, glycanJson, deisotoping);

            SearchMetaData searchInfo = new SearchMetaData();

            // targets
            while (tasks.TryDequeue(out SearchTask task))
            {
                TaskLocalSearch(ref tempResults, task,
                    precursorMatch, glycanSearch, searchInfo);

                searchCounter.Add(taskSize);
            }

            // decoys
            while (decoyTasks.TryDequeue(out SearchTask task))
            {
                TaskLocalSearch(ref tempDecoyResults, task,
                    precursorMatch, glycanSearch, searchInfo);

                searchCounter.Add(taskSize);
            }

            UpdateTask(tempResults, tempDecoyResults);
        }

    }
}
