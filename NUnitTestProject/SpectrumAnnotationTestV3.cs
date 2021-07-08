using MultiGlycanClassLibrary.util.mass;
using SpectrumProcess.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.annotation;
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
using SpectrumProcess.deisotoping;
using MultiGlycanTDLibrary.engine.score;

namespace NUnitTestProject
{
    public class SpectrumAnnotationTestV3
    {
        object obj = new object();
        string TypeToString(FragmentType type)
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
            string path = @"C:\Users\iruiz\Downloads\MSMS\134144_31_C18_120min_60oC_50cm.raw";
            string database = @"C:\Users\iruiz\Downloads\MSMS\database.json";
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
            ISearch<GlycanAnnotated> searcher3 = new BucketSearch<GlycanAnnotated>(
               ToleranceBy.Dalton, 0.5);
            GlycanAnnotationLazy annotator = new GlycanAnnotationLazy(searcher3,
                glycanJson.Parameters);
            CompdJson compdJson = glycanJson.Compound;

            // search
            double searchRange = 1.0;
            LocalNeighborPicking picking = new LocalNeighborPicking();
            IProcess process = new WeightedAveraging(new LocalNeighborPicking());
            Dictionary<int, List<PeakAnnotated>> final = new Dictionary<int, List<PeakAnnotated>>();

            ISearch<string> searcher = new BucketSearch<string>(ToleranceBy.PPM, 10);
            GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, compdJson);
            ISearch<Dictionary<FragmentType, List<string>>> searcher2 = 
                new BucketSearch<Dictionary<FragmentType, List<string>>>(ToleranceBy.Dalton, 0.5);
            
            Averagine averagine = new Averagine(AveragineType.PermethylatedGlycan);
            AveragineDeisotoping deisotoping = new AveragineDeisotoping(averagine,
                4, ToleranceBy.Dalton, 0.1);
            IGlycanSearch glycanSearch
                = new GlycanSearchDeisotoping(searcher2, glycanJson, deisotoping);
         
            SearchMetaData analyzer = new SearchMetaData();

            int targetScan = 18161;
            foreach (var scanPair in scanGroup)
            {
                if (scanPair.Value.Count > 0)
                {
                    int scan1 = scanPair.Key;
                    ISpectrum ms1 = reader.GetSpectrum(scan1);

                    foreach (int scan in scanPair.Value)
                    {
                        if (scan != targetScan)
                            continue;

                        double mz = reader.GetPrecursorMass(scan, reader.GetMSnOrder(scan));
                        List<IPeak> ms1Peaks = FilterPeaks(ms1.GetPeaks(), mz, searchRange);
                        if (ms1Peaks.Count() == 0)
                            continue;

                        ICharger charger = new Patterson();
                        int charge = charger.Charge(ms1Peaks, mz - searchRange, mz + searchRange);

                        // search
                        ISpectrum ms2 = reader.GetSpectrum(scan);
                        if (ms2.GetPeaks().Count <= 30)
                            continue;
                        ms2 = process.Process(ms2);

                        List<string> candidates = precursorMatch.Match(mz, charge);
                        if (candidates.Count == 0)
                            continue;
                        List<SearchResult> searched = glycanSearch.Search(candidates, ms2.GetPeaks(), charge, 1.0078);
                        List<SearchResult> results = analyzer.Commit(searched, mz, charge, scan, ms2.GetRetention());
                        List<PeakAnnotated> annotateds
                            = new List<PeakAnnotated>();

                        ClusterKMeans<IPeak>  cluster = new ClusterKMeans<IPeak>();

                        List<Point<IPeak>> points =
                            ms2.GetPeaks().Select(p => new Point<IPeak>(p.GetIntensity(), p)).ToList();
                        cluster.Run(points);
                        double minClusterIntensity = int.MaxValue;
                        int minClusterIndex = 0;
                        foreach (int index in cluster.Clusters.Keys)
                        {
                            double average =
                                cluster.Clusters[index].Average(peaks => peaks.Content().GetIntensity());
                            if (average < minClusterIntensity)
                            {
                                minClusterIntensity = average;
                                minClusterIndex = index;
                            }
                        }

                        foreach (SearchResult result in results)
                        {
                            int nTotal = result.Matches.Count;
                            int nMatched = 0;
                            foreach (int index in result.Matches.Keys)
                            {
                                int clusterIndex = cluster.Index[index];
                                if (clusterIndex == minClusterIndex)
                                    continue;
                                nMatched++;
                            }
                            result.Coverage = nMatched * 1.0 / nTotal;
                        }

                        foreach (SearchResult result in results)
                        {

                            result.Score = GlycanScorerHelper.ComputeScore(result, ms2.GetPeaks());
                            result.Fit = GlycanScorerHelper.ComputeFit(result, ms2.GetPeaks());
                        }

                        foreach (SearchResult result in results)
                        {
                            annotateds.AddRange(annotator.Annotated(ms2.GetPeaks(), result));
                        }


                        final[scan] = annotateds;
                    }

                }
            }



            //write out
            string outputPath = @"C:\Users\iruiz\Downloads\MSMS\annotated_spec.csv";
            //MultiGlycanClassLibrary.util.mass.Glycan.To.SetPermethylation(true, true);
            using (FileStream ostrm = new FileStream(outputPath, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan,mz,intensity,glycan,fragments");
                    string output = "";
                    foreach (var pair in final.OrderBy(p => p.Key))
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

    }
}
