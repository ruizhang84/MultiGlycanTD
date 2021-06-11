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
        List<SearchResult> decoys = new List<SearchResult>();
        private readonly object queueLock = new object();
        private readonly object resultLock = new object();
        private readonly double searchRange = 1;
        int taskSize = 0;
        int minPeaks = 30;

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
            tasks = new Queue<SearchTask>();
            GenerateTasks();
            taskSize = tasks.Count;
        }

        public List<SearchResult> Target()
            { return targets; }

        public List<SearchResult> Decoy()
            { return decoys; }

        SearchTask TryGetTask()
        {
            lock (queueLock)
            {
                if (tasks.Count > 0)
                    return tasks.Dequeue();
                return null;
            }
        }

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
                searches.Add(Task.Run(() => Search()));
            }

            Task.WaitAll(searches.ToArray());
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

                    //if (spectrum.Activation() != TypeOfMSActivation.CID)
                    //    continue;
                    SearchTask searchTask = new SearchTask(spectrum,
                        spectrum.PrecursorMZ(), spectrum.PrecursorCharge());
                    tasks.Enqueue(searchTask);
                }
            }

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
                                ISpectrum ms2 = reader.GetSpectrum(i);
                                if (ms2.GetPeaks().Count <= minPeaks)
                                    continue;
                                ms2 = process.Process(ms2);
                                SearchTask searchTask = new SearchTask(ms2, mz, charge);
                                tasks.Enqueue(searchTask);
                                
                            }
                        }
                            
                    }
                    readingCounter.Add(scanGroup.Count);
                });
            

        }

        void TaskSearch(ref List<SearchResult> results, 
            ref List<SearchResult> decoyResults,
            SearchTask task,
            GlycanPrecursorMatch precursorMatch,
            GlycanSearch glycanSearch,
            SearchAnalyzer searchAnalyzer)
        {
            foreach(double ion in SearchingParameters.Access.Ions)
            {
                //precursor match
                List<string> candidates = precursorMatch.Match(task.PrecursorMZ, task.Charge, ion);
                if (candidates.Count > 0)
                {
                    // spectrum search
                    List<SearchResult> searched = glycanSearch.Search(
                        task.Spectrum.GetPeaks(), task.Charge, candidates, false, ion);

                    if (searched.Count > 0)
                    {
                        // add meta data
                        List<SearchResult> temp = searchAnalyzer.Analyze(
                            searched, task.PrecursorMZ,
                            task.Spectrum.GetScanNum(),
                            task.Spectrum.GetRetention());
                        results.AddRange(temp);
                    }

                    // decoy spectrum search
                    List<SearchResult> decoySearched = glycanSearch.Search(
                        task.Spectrum.GetPeaks(), task.Charge, candidates, true, ion);

                    if (decoySearched.Count > 0)
                    {
                        // add meta data
                        List<SearchResult> temp = searchAnalyzer.Analyze(
                            decoySearched, task.PrecursorMZ,
                            task.Spectrum.GetScanNum(),
                            task.Spectrum.GetRetention());
                        decoyResults.AddRange(temp);
                    }
                }
            }
        }

        void Search()
        {
            List<SearchResult> tempResults = new List<SearchResult>();
            List<SearchResult> tempDecoyResults = new List<SearchResult>();

            ISearch<string> searcher = new BucketSearch<string>(
                SearchingParameters.Access.MS1ToleranceBy, 
                SearchingParameters.Access.MS1Tolerance);
            GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, compdJson, 
                SearchingParameters.Access.Cutoff);
          

            ISearch<string> searcher2 = new BucketSearch<string>(
                SearchingParameters.Access.MS2ToleranceBy, 
                SearchingParameters.Access.MSMSTolerance);
            GlycanSearch glycanSearch = new GlycanSearch(searcher2, glycanJson);

            SearchAnalyzer searchAnalyzer = new SearchAnalyzer();

            SearchTask task;
            while ((task = TryGetTask()) != null)
            {
                TaskSearch(ref tempResults, ref tempDecoyResults, task,
                    precursorMatch, glycanSearch, searchAnalyzer);

                searchCounter.Add(taskSize);
            }

            UpdateTask(tempResults, tempDecoyResults);
        }
    }
}
