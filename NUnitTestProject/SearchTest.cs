using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumProcess;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
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
            string path = @"C:\Users\Rui Zhang\Downloads\HBS1_dextrinspkd_C18_10252018.raw";
            string database = @"C:\Users\Rui Zhang\Downloads\small_database.json";
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
            string jsonStringRead = File.ReadAllText(database);
            GlycanJson glycanJson = JsonSerializer.Deserialize<GlycanJson>(jsonStringRead);
            CompdJson compdJson = glycanJson.Compound;


            // search
            double searchRange = 1.0;
            LocalNeighborPicking picking = new LocalNeighborPicking();
            IProcess process = new WeightedAveraging(new LocalNeighborPicking());

            object obj = new object();
            List<SearchResult> final = new List<SearchResult>();

            //Parallel.ForEach(scanGroup, scanPair =>
            foreach(var scanPair in scanGroup)
            {
                if (scanPair.Value.Count > 0)
                {
                    var watch = new System.Diagnostics.Stopwatch();
                    watch.Start();

                    int scan1 = scanPair.Key;
                    ISpectrum ms1 = reader.GetSpectrum(scan1);

                    foreach (int scan in scanPair.Value)
                    {
                        if (scan != 192)
                            continue;
                        double mz = reader.GetPrecursorMass(scan, reader.GetMSnOrder(scan));
                        List<IPeak> ms1Peaks = FilterPeaks(ms1.GetPeaks(), mz, searchRange);
                        ms1Peaks = picking.Process(ms1Peaks);


                        ICharger charger = new Patterson();
                        int charge = charger.Charge(ms1Peaks, mz - searchRange, mz + searchRange);

                        // search
                        ISpectrum ms2 = reader.GetSpectrum(scan);
                        ms2 = process.Process(ms2);

                        ISearch<string> searcher = new BucketSearch<string>(ToleranceBy.PPM, 10);
                        GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, compdJson, 0.01);
                        List<string> candidates = precursorMatch.Match(mz, charge);

                        ISearch<string> searcher2 = new BucketSearch<string>(ToleranceBy.Dalton, 0.01);
                        GlycanSearch glycanSearch = new GlycanSearch(searcher2, glycanJson);
                        List<SearchResult> searched = glycanSearch.Search(ms2.GetPeaks(), charge, candidates);

                        SearchAnalyzer analyzer = new SearchAnalyzer();
                        List<SearchResult> results = analyzer.Analyze(searched, mz, scan, ms2.GetRetention());

                        //EnvelopeProcess envelopeProcess = new EnvelopeProcess(ToleranceBy.Dalton, 0.01);
                        //GlycanEnvelopeMatch envelopeMatch = new GlycanEnvelopeMatch(envelopeProcess, compdJson);
                        //results = envelopeMatch.Match(results, ms1Peaks, mz, charge);


                        lock (obj)
                        {
                            final.AddRange(results);
                        }
                        //return;
                    }

                }
            }
        
            

            //write out
            string outputPath = @"C:\Users\Rui Zhang\Downloads\searching.csv";
            //MultiGlycanClassLibrary.util.mass.Glycan.To.SetPermethylation(true, true);
            using (FileStream ostrm = new FileStream(outputPath, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan,retention,glycan,mz,score,fit");
                    foreach(SearchResult r in final.OrderBy(p => p.Scan()))
                    {
                        string output = r.Scan().ToString() + ","
                            + r.Retention().ToString() + ","
                            + r.Glycan() + ","
                            + r.MZ().ToString() + ","
                            + r.Score().ToString() + ","
                            + r.Fit().ToString();
                        writer.WriteLine(output);
                    }
                    writer.Flush();
                }
            }

        }

    }
}
