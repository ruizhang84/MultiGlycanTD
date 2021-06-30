using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.annotation;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
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
using System.Text.Json;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumAnnotationTestV4
    {
        object obj = new object();

        List<FragmentTypes> types = new List<FragmentTypes>()
        {
            FragmentTypes.B, FragmentTypes.C, FragmentTypes.Y, FragmentTypes.Z,
            FragmentTypes.BY, FragmentTypes.BZ, FragmentTypes.CY, FragmentTypes.YY,
            FragmentTypes.YZ, FragmentTypes.ZZ
        };

        string TypeToString(FragmentTypes type)
        {
            switch (type)
            {
                case FragmentTypes.B:
                    return "B";
                case FragmentTypes.C:
                    return "C";
                case FragmentTypes.Y:
                    return "Y";
                case FragmentTypes.Z:
                    return "Z";
                case FragmentTypes.BY:
                    return "BY";
                case FragmentTypes.BZ:
                    return "BZ";
                case FragmentTypes.CY:
                    return "CY";
                case FragmentTypes.YY:
                    return "YY";
                case FragmentTypes.YZ:
                    return "YZ";
                case FragmentTypes.ZZ:
                    return "ZZ";
            }
            return "";
        }

        List<IGlycan> FragmentsBuild(FragmentTypes type, IGlycan glycan)
        {
            switch (type)
            {
                case FragmentTypes.B:
                case FragmentTypes.C:
                    return GlycanFragmentBuilder.BionsLikeFragments(glycan);
                case FragmentTypes.Y:
                case FragmentTypes.Z:
                    return GlycanFragmentBuilder.YionsLikeFragments(glycan);
                case FragmentTypes.BY:
                case FragmentTypes.BZ:
                case FragmentTypes.CY:
                    return GlycanFragmentBuilder.BYionsLikeFragments(glycan);
                case FragmentTypes.YY:
                case FragmentTypes.YZ:
                case FragmentTypes.ZZ:
                    return GlycanFragmentBuilder.YYionsLikeFragments(glycan);

            }
            return new List<IGlycan>();
        }

        double MassBuild(FragmentTypes type, IGlycan glycan)
        {
            switch (type)
            {
                case FragmentTypes.B:
                    return GlycanIonsBuilder.Build.Bion(glycan);
                case FragmentTypes.C:
                    return GlycanIonsBuilder.Build.Cion(glycan);
                case FragmentTypes.Y:
                    return GlycanIonsBuilder.Build.Yion(glycan);
                case FragmentTypes.Z:
                    return GlycanIonsBuilder.Build.Zion(glycan);
                case FragmentTypes.BY:
                    return GlycanIonsBuilder.Build.BYion(glycan);
                case FragmentTypes.BZ:
                    return GlycanIonsBuilder.Build.BZion(glycan);
                case FragmentTypes.CY:
                    return GlycanIonsBuilder.Build.CYion(glycan);
                case FragmentTypes.YY:
                    return GlycanIonsBuilder.Build.YYion(glycan);
                case FragmentTypes.YZ:
                    return GlycanIonsBuilder.Build.YZion(glycan);
                case FragmentTypes.ZZ:
                    return GlycanIonsBuilder.Build.ZZion(glycan);
            }
            return 0;
        }

        void BuildMassMap(string id, IGlycan glycan, FragmentTypes type,
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
            string path = @"D:\Data\peptide_standard.mgf";
            string database = @"C:\Users\iruiz\Downloads\MSMS\database.json";
            MGFSpectrumReader mgfReader = new MGFSpectrumReader();
            mgfReader.Init(path);


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
                foreach (FragmentTypes type in types)
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
            ISearch<Dictionary<FragmentTypes, List<string>>> searcher2 = new BucketSearch<Dictionary<FragmentTypes, List<string>>>(ToleranceBy.Dalton, 0.1);
            GlycanSearch glycanSearch = new GlycanSearch(searcher2, glycanJson);
            ISearch<GlycanAnnotated> searcher3 = new BucketSearch<GlycanAnnotated>(ToleranceBy.Dalton, 0.1);
            SearchMetaData analyzer = new SearchMetaData();
            GlycanAnnotation glycanAnnotation = new GlycanAnnotation(searcher3,
                massMap.ToDictionary(entry => entry.Key, entry => entry.Value));

            int targetScan = 55488;
            double targetMZ = 447.9912;
            double delta = 0; //  809.428345 - 799.423218;


            Dictionary<int, MS2Spectrum> spectraData = mgfReader.GetSpectrum();
            foreach (int scan in spectraData.Keys)
            {
                if (scan != targetScan)
                    continue;

                MS2Spectrum ms2 = spectraData[scan];

                if (ms2.GetPeaks().Count <= 30)
                    continue;
                ms2 = process.Process(ms2) as MS2Spectrum;
                foreach (IPeak pk in ms2.GetPeaks())
                {
                    pk.SetMZ(pk.GetMZ() + delta);
                }
                double mz = ms2.PrecursorMZ();
                int charge = ms2.PrecursorCharge();

                if (targetMZ > 0)
                    mz = targetMZ;
                List<string> candidates = precursorMatch.Match(mz, charge);
                if (candidates.Count == 0)
                    continue;
                List<SearchResult> searched = glycanSearch.Search(ms2.GetPeaks(), charge, candidates);
                List<SearchResult> results = analyzer.Commit(searched, mz, charge, scan, ms2.GetRetention());
                List<PeakAnnotated> annotateds = glycanAnnotation.Annotated(ms2.GetPeaks(), charge, results);

                final[scan] = annotateds;
            }

            //write out
            string outputPath = @"C:\Users\iruiz\Downloads\MSMS\annotated_spec2"
                    + (targetMZ > 0 ? "_decoy" : "_target") + ".csv";
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
