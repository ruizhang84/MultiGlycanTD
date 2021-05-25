using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumData.Spectrum;
using SpectrumProcess;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTD
{

    public class MultiThreadingSearch
    {
        Counter readingCounter;
        Counter searchCounter;
        GlycanJson glycanJson;
        CompdJson compdJson;

        Queue<SearchTask> tasks;
        string msPath;
        List<SearchResult> targets = new List<SearchResult>();
        private readonly object queueLock = new object();
        private readonly object resultLock = new object();
        private readonly double searchRange = 1;
        int taskSize = 0;
        double cutoff = 0.01;

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
            tasks = new Queue<SearchTask>();
            GenerateTasks();
            taskSize = tasks.Count;
        }

        public List<SearchResult> Target()
            { return targets; }

        SearchTask TryGetTask()
        {
            lock (queueLock)
            {
                if (tasks.Count > 0)
                    return tasks.Dequeue();
                return null;
            }
        }

        void UpdateTask(List<SearchResult> t)
        {
            lock (resultLock)
            {
                targets.AddRange(t);
            }
        }

        public void Run()
        {
            List<Task> searches = new List<Task>();
            for (int i = 0; i < SearchingParameters.Access.ThreadNums; i++)
            {
                searches.Add(Task.Run(() => Search()));
            }

            Task.WaitAll(searches.ToArray());
        }

        void GenerateTasks()
        {
            // read spectrum
            ISpectrumReader reader = new ThermoRawSpectrumReader();
            IProcess picking = new LocalNeighborPicking();
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
                                ms1Peaks = picking.Process(ms1Peaks);

                                // charage
                                ICharger charger = new Patterson();
                                int charge = charger.Charge(ms1Peaks, mz - searchRange, mz + searchRange);


                                // search
                                ISpectrum ms2 = reader.GetSpectrum(scan);
                                ms2 = process.Process(ms2);

                                if (ms2.GetPeaks().Count > 0)
                                {
                                    SearchTask searchTask = new SearchTask(ms2, ms1Peaks, mz, charge);
                                    tasks.Enqueue(searchTask);
                                }
                            }
                        }
                            
                    }
                    readingCounter.Add(scanGroup.Count);
                });
            

        }

        void TaskSearch(ref List<SearchResult> results, SearchTask task,
            GlycanPrecursorMatch precursorMatch,
            GlycanEnvelopeMatch envelopeMatch,
            GlycanSearch glycanSearch,
            SearchAnalyzer searchAnalyzer)
        {
            // get spectrum
            ISpectrum spectrum = task.Spectrum;
            //precursor match
            List<string> candidates = precursorMatch.Match(task.PrecursorMZ, task.Charge);
            if (candidates.Count > 0)
            {
                // spectrum search
                List<SearchResult> searched = glycanSearch.Search(
                    task.Spectrum.GetPeaks(), task.Charge, candidates);

                // isotopic cluster match
                List<SearchResult> matched = envelopeMatch.Match(searched,
                    task.MSPeaks, task.PrecursorMZ, task.Charge);

                // add meta data
                List < SearchResult> temp = searchAnalyzer.Analyze(
                    matched, task.PrecursorMZ,
                    task.Spectrum.GetScanNum(),
                    task.Spectrum.GetRetention());
                results.AddRange(temp);
            }
        }

        void Search()
        {
            List<SearchResult> tempResults = new List<SearchResult>();


            ISearch<string> searcher = new BucketSearch<string>(
                SearchingParameters.Access.MS1ToleranceBy, 
                SearchingParameters.Access.MS1Tolerance);
            GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, compdJson, cutoff);
          

            ISearch<string> searcher2 = new BucketSearch<string>(
                SearchingParameters.Access.MS2ToleranceBy, 
                SearchingParameters.Access.MSMSTolerance);
            GlycanSearch glycanSearch = new GlycanSearch(searcher2, glycanJson);

            SearchAnalyzer searchAnalyzer = new SearchAnalyzer();

            EnvelopeProcess envelopeProcess = new EnvelopeProcess(
                SearchingParameters.Access.MS2ToleranceBy,
                SearchingParameters.Access.MSMSTolerance);
            GlycanEnvelopeMatch envelopeMatch = new GlycanEnvelopeMatch(envelopeProcess, compdJson);

            SearchTask task;
            while ((task = TryGetTask()) != null)
            {
                TaskSearch(ref tempResults, task,
                    precursorMatch, envelopeMatch,
                    glycanSearch, searchAnalyzer);

                searchCounter.Add(taskSize);
            }

            UpdateTask(tempResults);
        }
    }
}
