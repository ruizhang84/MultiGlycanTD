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

namespace NUnitTestProject
{
    public class SpectrumAnnotationTestV3
    {
        object obj = new object();

        List<FragmentType> types = new List<FragmentType>()
        {
            FragmentType.B, FragmentType.Y,
            FragmentType.BY, FragmentType.YY
        };

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

        List<IGlycan> FragmentsBuild(FragmentType type, IGlycan glycan)
        {
            switch (type)
            {
                case FragmentType.B:
                case FragmentType.C:
                    return GlycanFragmentBuilder.BionsLikeFragments(glycan);
                case FragmentType.Y:
                case FragmentType.Z:
                    return GlycanFragmentBuilder.YionsLikeFragments(glycan);
                case FragmentType.BY:
                case FragmentType.BZ:
                case FragmentType.CY:
                    return GlycanFragmentBuilder.BYionsLikeFragments(glycan);
                case FragmentType.YY:
                case FragmentType.YZ:
                case FragmentType.ZZ:
                    return GlycanFragmentBuilder.YYionsLikeFragments(glycan);

            }
            return new List<IGlycan>();
        }

        double MassBuild(FragmentType type, IGlycan glycan)
        {
            switch (type)
            {
                case FragmentType.B:
                    return GlycanIonsBuilder.Build.Bion(glycan);
                case FragmentType.C:
                    return GlycanIonsBuilder.Build.Cion(glycan);
                case FragmentType.Y:
                    return GlycanIonsBuilder.Build.Yion(glycan);
                case FragmentType.Z:
                    return GlycanIonsBuilder.Build.Zion(glycan);
                case FragmentType.BY:
                    return GlycanIonsBuilder.Build.BYion(glycan);
                case FragmentType.BZ:
                    return GlycanIonsBuilder.Build.BZion(glycan);
                case FragmentType.CY:
                    return GlycanIonsBuilder.Build.CYion(glycan);
                case FragmentType.YY:
                    return GlycanIonsBuilder.Build.YYion(glycan);
                case FragmentType.YZ:
                    return GlycanIonsBuilder.Build.YZion(glycan);
                case FragmentType.ZZ:
                    return GlycanIonsBuilder.Build.ZZion(glycan);
            }
            return 0;
        }

        void BuildMassMap(string id, IGlycan glycan, FragmentType type,
            ref ConcurrentDictionary<double, List<GlycanAnnotated>> massMap)
        {
            if (glycan.IsValid())
            {
                List<IGlycan> bionsLikeFragments = FragmentsBuild(type, glycan);
                foreach (IGlycan g in bionsLikeFragments)
                {
                    double mass = MassBuild(type, g);
                    lock (obj)
                    {
                        if (!massMap.ContainsKey(mass))
                        {
                            massMap[mass] = new List<GlycanAnnotated>();
                        }
                        massMap[mass].Add(new GlycanAnnotated
                        {
                            Parent = id,
                            Type = type,
                            Glycan = g.ID()
                        });
                    }
                }
            }
        }

        [Test]
        public void SearchSpectrum()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
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
            CompdJson compdJson = glycanJson.Compound;

            GlycanBuilder glycanBuilder = new GlycanBuilder(7, 7, 5, 4, 0,
                true, true, true,
                10, true, true);
            glycanBuilder.Build();
            var map = glycanBuilder.GlycanMaps();

            GlycanIonsBuilder.Build.Permethylated = true;
            GlycanIonsBuilder.Build.Reduced = true;
            Glycan.To.SetPermethylation(true, true);

            ConcurrentDictionary<double, List<GlycanAnnotated>> massMap
                = new ConcurrentDictionary<double, List<GlycanAnnotated>>();
            Parallel.ForEach(map, pair =>
            {
                foreach (FragmentType type in types)
                {
                    BuildMassMap(pair.Key, pair.Value, type, ref massMap);
                }
            });

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


            ISearch<GlycanAnnotated> searcher3 = 
                new BucketSearch<GlycanAnnotated>(ToleranceBy.Dalton, 0.5);
            SearchMetaData analyzer = new SearchMetaData();
            GlycanAnnotation glycanAnnotation = new GlycanAnnotation(searcher3,
                massMap.ToDictionary(entry => entry.Key, entry => entry.Value));

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
                        ISpectrum ms2 = reader.GetSpectrum(scan);
                        if (ms2.GetPeaks().Count <= 30)
                            continue;
                        ms2 = process.Process(ms2);

                        List<string> candidates = precursorMatch.Match(mz, charge);
                        if (candidates.Count == 0)
                            continue;
                        List<SearchResult> searched = glycanSearch.Search(candidates, ms2.GetPeaks(), charge, 1.0078);
                        List<SearchResult> results = analyzer.Commit(searched, mz, charge, scan, ms2.GetRetention());

                        List<PeakAnnotated> annotateds = glycanAnnotation.Annotated(deisotoping.Process(ms2.GetPeaks(), 1.0078), charge, results);

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
