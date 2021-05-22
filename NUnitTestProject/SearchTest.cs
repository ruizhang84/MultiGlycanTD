using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.engine.search;
using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumProcess;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SearchTest
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
        public void SearchSpectrum()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\HBS1_dextrinspkd_C18_10252018.raw";
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

            // init
            GlycanBuilder glycanBuilder = new GlycanBuilder();
            glycanBuilder.Build();


            // search
            double searchRange = 1.0;
            IProcess picking = new LocalNeighborPicking();
            IProcess process = new WeightedAveraging(new LocalNeighborPicking());
            foreach (var scanPair in scanGroup)
            {
                if (scanPair.Value.Count > 0)
                {
                    int scan1 = scanPair.Key;
                    ISpectrum ms1 = reader.GetSpectrum(scanPair.Key);
                    ms1 = picking.Process(ms1);
                    List<IPeak> ms1Peaks = ms1.GetPeaks();

                    foreach (int scan in scanPair.Value)
                    {
                        double mz = reader.GetPrecursorMass(scan1, reader.GetMSnOrder(scan1));
                        if (FilterPeaks(ms1Peaks, mz, searchRange).Count == 0) 
                            continue;


                        ICharger charger = new Patterson();
                        int charge = charger.Charge(ms1.GetPeaks(), mz - searchRange, mz + searchRange);


                        // search
                        ISpectrum ms2 = reader.GetSpectrum(scan);
                        ms2 = process.Process(ms2);

                        ISearch<string> searcher = new BucketSearch<string>(ToleranceBy.PPM, 10);
                        GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, 
                            glycanBuilder.GlycanCompositionMaps(), glycanBuilder.GlycanMasMaps());

                    }





                }
            }

        }

    }
}
