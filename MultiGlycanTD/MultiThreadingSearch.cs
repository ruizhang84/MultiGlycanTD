using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
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

namespace MultiGlycanTD
{

    public class MultiThreadingSearch
    {
        Counter readingCounter;
        Counter searchCounter;
        GlycanJson glycanJson;
        CompdJson compdJson;

        ConcurrentQueue<SearchTask> tasks;
        ConcurrentQueue<SearchTask> decoyTasks;
        string msPath;
        List<SearchResult> targets = new List<SearchResult>();
        List<SearchResult> decoys = new List<SearchResult>();
        private readonly object resultLock = new object();
        private readonly double searchRange = 1;
        int taskSize = 0;
        // It is not likely that a target spectrum has very high charge
        // for glycan or very few peaks for fragments.
        int minPeaks = 30;   // sequencable spectrum with min num peaks
        int maxCharge = 3;   // max charge of spectrum to consider
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
                Task LastTask = new Task(() => Search());
                LastTask.Start();
                searches.Add(LastTask);
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
                                if (charge > maxCharge)
                                    continue;

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

        void GenerateDecoyTasks()
        {
            // find max precursor charges
            int maxCharge = 0;
            foreach (SearchTask task in tasks)
            {
                maxCharge = Math.Max(maxCharge, task.Charge);
            }

            // split spectrum on charges
            Parallel.For(1, maxCharge + 1, charge =>
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
              });
        }

        void TaskSearch(ref List<SearchResult> results,
            SearchTask task,
            GlycanPrecursorMatch precursorMatch,
            GlycanSearch glycanSearch,
            SearchMetaData searchInfo)
        {
            foreach (double ion in SearchingParameters.Access.Ions)
            {
                //precursor match
                List<string> candidates = precursorMatch.Match(task.PrecursorMZ, task.Charge, ion);
                if (candidates.Count > 0)
                {
                    // spectrum search
                    List<SearchResult> searched = glycanSearch.Search(
                        task.Spectrum.GetPeaks(), task.Charge, candidates, ion);

                    if (searched.Count > 0)
                    {
                        // add meta data
                        List<SearchResult> temp = searchInfo.Commit(
                            searched, task.PrecursorMZ,
                            task.Spectrum.GetScanNum(),
                            task.Spectrum.GetRetention());
                        results.AddRange(temp);
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

            SearchMetaData searchInfo = new SearchMetaData();

            // targets
            while (tasks.TryDequeue(out SearchTask task))
            {
                TaskSearch(ref tempResults, task,
                    precursorMatch, glycanSearch, searchInfo);

                searchCounter.Add(taskSize);
            }

            // decoys
            while (decoyTasks.TryDequeue(out SearchTask task))
            {
                TaskSearch(ref tempDecoyResults, task,
                    precursorMatch, glycanSearch, searchInfo);

                searchCounter.Add(taskSize);
            }

            UpdateTask(tempResults, tempDecoyResults);
        }

    }
}
