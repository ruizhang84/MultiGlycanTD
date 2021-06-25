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
    public class SpectrumSimUnitTest
    {
        public double computeDot(List<IPeak> l1, List<IPeak> l2)
        {
            if (l1.Count == 0 || l2.Count == 0)
            {
                return 0;
            }
            int count = Math.Min(l1.Count, l2.Count);
            List<IPeak> t1 = l1.OrderByDescending(x => x.GetIntensity())
                .Take(count).OrderBy(x => x.GetMZ()).ToList();
            List<IPeak> t2 = l2.OrderByDescending(x => x.GetIntensity()).
                Take(count).OrderBy(x => x.GetMZ()).ToList();
            double numerator = 0;
            for (int i = 0; i < t1.Count; i++)
            {
                numerator += t1[i].GetIntensity() * t2[i].GetIntensity();
            }
            return numerator;
        }

        public double computeCos(List<IPeak> p1, List<IPeak> p2, double tol)
        {
            double lowerBound = Math.Min(p1.Min(x => x.GetMZ()), p2.Min(x => x.GetMZ()));
            double upperBound = Math.Max(p1.Max(x => x.GetMZ()), p2.Max(x => x.GetMZ()));
            int bucketNums = (int)Math.Ceiling((upperBound - lowerBound + 1) / tol);

            List<IPeak>[] q1 = new List<IPeak>[bucketNums];
            List<IPeak>[] q2 = new List<IPeak>[bucketNums];

            for (int i = 0; i < bucketNums; i++)
            {
                q1[i] = new List<IPeak>();
                q2[i] = new List<IPeak>();
            }

            foreach (IPeak pk in p1)
            {
                int index = (int)Math.Ceiling((pk.GetMZ() - lowerBound) / tol);
                q1[index].Add(pk);
            }
            foreach (IPeak pk in p2)
            {
                int index = (int)Math.Ceiling((pk.GetMZ() - lowerBound) / tol);
                q2[index].Add(pk);
            }

            double numerator = 0;
            for (int i = 0; i < bucketNums; i++)
            {
                numerator += computeDot(q1[i], q2[i]);
            }

            double denominator1 = 0;
            foreach (IPeak pk in p1)
            {
                denominator1 += pk.GetIntensity() * pk.GetIntensity();
            }
            double denominator2 = 0;
            foreach (IPeak pk in p2)
            {
                denominator2 += pk.GetIntensity() * pk.GetIntensity();
            }
            double denominator = Math.Sqrt(denominator1) * Math.Sqrt(denominator2);
            return numerator / denominator;
        }
        [Test]
        public void CosTest()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
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
            GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, compdJson, 0.01);

            List<int> scanMSMS = new List<int>();
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
                        scanMSMS.Add(scan);
                    }

                }
            }

            string outputPath = @"C:\Users\iruiz\Downloads\MSMS\sim.csv";
            using (FileStream ostrm = new FileStream(outputPath, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    Dictionary<int, ISpectrum> mem = new Dictionary<int, ISpectrum>();
                    for (int i = 0; i < scanMSMS.Count; i++)
                    {
                        ISpectrum A = reader.GetSpectrum(scanMSMS[i]);
                        if (A.GetPeaks().Count > 0)
                            mem[i] = A;
                    }

                    string output = "";
                    scanMSMS = mem.Keys.ToList();
                    writer.WriteLine("scan1,scan2,cos");
                    for (int i = 0; i < scanMSMS.Count - 1; i++)
                    {
                        for (int j = i + 1; j < scanMSMS.Count; j++)
                        {
                            ISpectrum A = mem[scanMSMS[i]];
                            ISpectrum B = mem[scanMSMS[j]];
                            double sim = computeCos(A.GetPeaks(), B.GetPeaks(), 0.1);
                            output += scanMSMS[i].ToString() + "," + scanMSMS[j].ToString() + "," + sim.ToString() + "\n";

                        }
                    }
                    writer.WriteLine(output);
                }
            }
                   
            Assert.Pass();
        }

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
    }
}
