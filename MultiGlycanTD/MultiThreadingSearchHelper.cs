using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.search;
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
    public class ProgressingEventArgs : EventArgs
    {
        public int Total { get; set; }
    }

    public class Counter
    {
        public event EventHandler<ProgressingEventArgs> progressChange;

        protected virtual void OnProgressChanged(ProgressingEventArgs e)
        {
            EventHandler<ProgressingEventArgs> handler = progressChange;
            handler?.Invoke(this, e);
        }

        public void Add(int total)
        {
            ProgressingEventArgs e = new ProgressingEventArgs
            {
                Total = total
            };
            OnProgressChanged(e);
        }
    }

    public class SearchTask
    {
        public SearchTask(ISpectrum spectrum, double mz, int charge)
        {
            Spectrum = spectrum;
            PrecursorMZ = mz;
            Charge = charge;
        }
        public ISpectrum Spectrum { get; set; }
        public double PrecursorMZ { get; set; }
        public int Charge { get; set; }

    }

    public class MultiThreadingSearchHelper
    {
        public static void GenerateSearchTasks(string spectrumPath,
           ConcurrentQueue<SearchTask> searchTasks,
           ConcurrentDictionary<int, ISpectrum> tandemSpectra,
           Counter readingCounter, 
           int minPeaks, int maxCharge, int minCharge, double searchRange)
        {
            if (Path.GetExtension(spectrumPath) == ".mgf")
            {
                MGFSpectrumReader mgfReader = new MGFSpectrumReader();
                mgfReader.Init(spectrumPath);

                Dictionary<int, MS2Spectrum> spectraData = mgfReader.GetSpectrum();
                foreach (int scan in spectraData.Keys)
                {
                    MS2Spectrum spectrum = spectraData[scan];
                    readingCounter.Add(spectraData.Count);
                    if (spectrum.GetPeaks().Count <= minPeaks)
                        continue;

                    tandemSpectra[scan] = spectrum;
                    SearchTask searchTask = new SearchTask(spectrum,
                        spectrum.PrecursorMZ(), spectrum.PrecursorCharge());
                    searchTasks.Enqueue(searchTask);
                }
            }

            // read spectrum
            ISpectrumReader reader = new ThermoRawSpectrumReader();
            IProcess process = new WeightedAveraging(new LocalNeighborPicking());

            reader.Init(spectrumPath);

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
                                if (charge > maxCharge || charge < minCharge)
                                    continue;

                                // search
                                ISpectrum ms2 = reader.GetSpectrum(i);
                                if (ms2.GetPeaks().Count <= minPeaks)
                                    continue;
                                ms2 = process.Process(ms2);
                                tandemSpectra[i] = ms2;
                                SearchTask searchTask = new SearchTask(ms2, mz, charge);
                                searchTasks.Enqueue(searchTask);

                            }
                        }

                    }
                    readingCounter.Add(scanGroup.Count);
                });
        }

        public static List<IPeak> FilterPeaks(List<IPeak> peaks, double target, double range)
        {
            if (peaks.Count == 0)
            {
                return peaks;
            }

            int start = 0;
            int end = peaks.Count - 1;
            int middle = 0;
            if (peaks[start].GetMZ() > target - range)
            {
                middle = start;
            }
            else
            {
                while (start + 1 < end)
                {
                    middle = (end - start) / 2 + start;
                    double mz = peaks[middle].GetMZ() + range;
                    if (mz == target)
                    {
                        break;
                    }
                    else if (mz < target)
                    {
                        start = middle;
                    }
                    else
                    {
                        end = middle - 1;
                    }
                }
            }

            List<IPeak> res = new List<IPeak>();
            while (middle < peaks.Count)
            {
                if (peaks[middle].GetMZ() > target + range)
                    break;
                res.Add(peaks[middle++]);
            }
            return res;
        }
        public static void Report(
            string path, 
            List<SearchResult> results)
        {
            using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan,retention,glycan,struct,precursor_mz,score");
                    foreach (SearchResult r in results)
                    {
                        string output = r.Scan.ToString() + ","
                            + r.Retention.ToString() + ","
                            + r.Composition + ","
                            + r.Glycan + ","
                            + r.PrecursorMZ.ToString() + ","
                            + r.Score.ToString();
                        writer.WriteLine(output);
                    }
                    writer.Flush();
                }
            }
        }

    }
}
