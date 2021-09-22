using MultiGlycanTDLibrary.engine.annotation;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumData.Spectrum;
using SpectrumProcess;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
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

        public SearchTask(ISpectrum spectrum, double mz, int charge, List<IPeak> peaks)
        {
            Spectrum = spectrum;
            PrecursorMZ = mz;
            Charge = charge;
            Peaks = peaks;
        }

        public ISpectrum Spectrum { get; set; }
        public double PrecursorMZ { get; set; }
        public int Charge { get; set; }
        public List<IPeak> Peaks { get; set; } = null; // ms1 peaks

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
                    if (spectrum.PrecursorCharge() > maxCharge)
                        continue;

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
                                tandemSpectra[i] = new MS2Spectrum(ms2, mz, charge);
                                SearchTask searchTask = new SearchTask(ms2, mz, charge, process.Process(ms1Peaks));
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
                            + r.Score.ToString()
                            ;
                        writer.WriteLine(output);
                    }
                    writer.Flush();
                }
            }
        }

        static string TypeToString(FragmentType type)
        {
            switch (type)
            {
                case FragmentType.B:
                    return "B";
                case FragmentType.C:
                    return "C";
                case FragmentType.Y:
                    return "Y";
                case FragmentType.Z:
                    return "Z";
                case FragmentType.BY:
                    return "BY";
                case FragmentType.BZ:
                    return "BZ";
                case FragmentType.CY:
                    return "CY";
                case FragmentType.YY:
                    return "YY";
                case FragmentType.YZ:
                    return "YZ";
                case FragmentType.ZZ:
                    return "ZZ";
                case FragmentType.BYY:
                    return "BYY";
                case FragmentType.BYZ:
                    return "BYZ";
                case FragmentType.BZZ:
                    return "BZZ";
                case FragmentType.CYY:
                    return "CYY";
                case FragmentType.CYZ:
                    return "CYZ";
                case FragmentType.CZZ:
                    return "CZZ";
                case FragmentType.YYY:
                    return "YYY";
                case FragmentType.YYZ:
                    return "YYZ";
                case FragmentType.YZZ:
                    return "YZZ";
                case FragmentType.ZZZ:
                    return "ZZZ";
            }
            return "";
        }

        public static void AnnotationReport(
            string path, Dictionary<int, List<PeakAnnotated>> annotation)
        {
            using (FileStream ostrm = new FileStream(path, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan,mz,intensity,glycan,fragments");
                    string output = "";
                    foreach (var pair in annotation.OrderBy(p => p.Key))
                    {
                        int scan = pair.Key;
                        List<PeakAnnotated> peakAnnotateds = pair.Value;
                        foreach (var pka in peakAnnotateds)
                        {
                            output += scan.ToString() + "," +
                                pka.Peak.GetMZ() + "," +
                                pka.Peak.GetIntensity() + "," +
                                pka.Glycan + "," +
                                string.Join("|", pka.Fragments.Select(f => TypeToString(f.Type) + ":" + f.Glycan)) + "\n";

                        }
                    }
                    writer.WriteLine(output);
                    writer.Flush();
                }
            }
        }

        public static Dictionary<string, List<double>> ReadGlycanDiagnosticPeaks(string path)
        {
            Dictionary<string, List<double>> glycanDiagnosticPeaks = 
                new Dictionary<string, List<double>>();

            using (FileStream fileStream = new FileStream(path, FileMode.Open, FileAccess.Read))
            {
                using (StreamReader sr = new StreamReader(fileStream))
                {
                    string line;

                    // Read lines from the file until end of file (EOD) is reached.
                    while ((line = sr.ReadLine()) != null)
                    {
                        string[] values = line.Split(',');
                        if (values.Length == 3)
                        {
                            string glycan = values[1];
                            if (double.TryParse(values[2], out double mz))
                            {
                                if (!glycanDiagnosticPeaks.ContainsKey(glycan))
                                {
                                    glycanDiagnosticPeaks[glycan] = new List<double>();
                                }
                                glycanDiagnosticPeaks[glycan].Add(mz);
                            }
                        }

                    }
                }
            }

            return glycanDiagnosticPeaks;
        }

    }
}
