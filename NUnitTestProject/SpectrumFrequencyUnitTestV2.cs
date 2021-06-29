using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumProcess;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumFrequencyUnitTestV2
    {
        List<IPeak> FilterPeaks(List<IPeak> peaks, double target, double range)
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

        [Test]
        public void CountPeaksTest()
        {
            // read spectrum
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
            string output = @"C:\Users\iruiz\Downloads\MSMS\peaks_glycan.csv";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

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

            string database = @"C:\Users\iruiz\Downloads\MSMS\database.json";
            string jsonStringRead = File.ReadAllText(database);
            GlycanJson glycanJson = JsonSerializer.Deserialize<GlycanJson>(jsonStringRead);
            CompdJson compdJson = glycanJson.Compound;

            double searchRange = 1.0;

            ISearch<string> searcher = new BucketSearch<string>(ToleranceBy.PPM, 10);
            GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, compdJson);

            Dictionary<double, int> counts = new Dictionary<double, int>();
            foreach (var scanPair in scanGroup)
            {
                if (scanPair.Value.Count > 0)
                {
                    int scan1 = scanPair.Key;
                    ISpectrum ms1 = reader.GetSpectrum(scan1);

                    foreach (int scan in scanPair.Value)
                    {

                        double mz = reader.GetPrecursorMass(scan, reader.GetMSnOrder(scan));
                        List<IPeak> ms1Peaks = FilterPeaks(ms1.GetPeaks(), mz, searchRange);
                        if (ms1Peaks.Count() == 0)
                            continue;

                        ICharger charger = new Patterson();
                        int charge = charger.Charge(ms1Peaks, mz - searchRange, mz + searchRange);

                        // search
                        List<string> candidates = precursorMatch.Match(mz, charge);
                        if (candidates.Count == 0)
                            continue;
                        foreach (IPeak peak in reader.GetSpectrum(scan).GetPeaks())
                        {
                            double mz2 = Math.Round(peak.GetMZ(), 1);
                            if (!counts.ContainsKey(mz2))
                            {
                                counts[mz2] = 0;
                            }
                            counts[mz2] += 1;
                        }
                    }

                }
            }

            using (FileStream ostrm = new FileStream(output, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("mz,#peaks");
                    foreach (double mz in counts.Keys)
                    {
                        writer.WriteLine(mz.ToString() + "," +
                            counts[mz].ToString());
                    }
                }
            }
        }

       
    }
}
