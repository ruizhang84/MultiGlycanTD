using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using MultiGlycanTDLibrary.model.glycan;
using NUnit.Framework;
using SpectrumData.Reader;
using SpectrumData.Spectrum;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumGlycanSearchUnitTest
    {
        [Test]
        public void SearchTest()
        {
            // init database
            GlycanBuilder glycanBuilder =
                new GlycanBuilder(7, 7, 5, 4, 0, true, false, false, 10, 
                false, false, Derivatization.Underivatized);
            glycanBuilder.Build();

            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();
            CompdJson compdJson = new CompdJson()
            {
                DistrMap = distr_map,
                MassMap = mass_map
            };

            object obj = new object();
            var map = glycanBuilder.GlycanMaps();
            Dictionary<double, List<string>> fragments = new Dictionary<double, List<string>>();
            List<Tuple<string, List<double>>> fragmentsContainer =
                new List<Tuple<string, List<double>>>();

            GlycanIonsBuilder.Build.Permethylated = false;
            GlycanIonsBuilder.Build.Derivatization = Glycan.kWater;
            Glycan.To.SetPermethylation(false, false);
            Glycan.To.Derivatization = Glycan.kWater;
            

            Parallel.ForEach(map, pair =>
            {
                var id = pair.Key;
                var glycan = pair.Value;
                if (glycan.IsValid())
                {
                    List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                                        .OrderBy(m => m).Select(m => Math.Round(m, 4)).ToList();
                    lock (obj)
                    {
                        fragmentsContainer.Add(Tuple.Create(id, massList));
                    }
                }
            });

            foreach (Tuple<string, List<double>> item in fragmentsContainer)
            {
                string id = item.Item1;
                foreach (double mass in item.Item2)
                {
                    if (!fragments.ContainsKey(mass))
                        fragments[mass] = new List<string>();
                    fragments[mass].Add(id);
                }
            }

            var composition_map = glycanBuilder.GlycanCompositionMaps();
            Dictionary<string, List<string>> id_map = new Dictionary<string, List<string>>();
            foreach (var pair in composition_map)
            {
                id_map[pair.Key] = pair.Value.Select(p => p.ID()).ToList();
            }

            GlycanJson glycanJson = new GlycanJson()
            {
                Compound = compdJson,
                IDMap = id_map,
                Fragments = fragments
            };

            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\glycan_standard.mgf";
  
            MGFSpectrumReader mgfReader = new MGFSpectrumReader();
            mgfReader.Init(path);


            // search
            List<SearchResult> final = new List<SearchResult>();

            ISearch<string> searcher = new BucketSearch<string>(ToleranceBy.Dalton, 1.0);
            GlycanPrecursorMatch precursorMatch = new GlycanPrecursorMatch(searcher, compdJson, 0.01);
            ISearch<string> searcher2 = new BucketSearch<string>(ToleranceBy.Dalton, 0.1);
            GlycanSearch glycanSearch = new GlycanSearch
                (searcher2, glycanJson);
            SearchAnalyzer analyzer = new SearchAnalyzer();

            Dictionary<int, MS2Spectrum> spectraData = mgfReader.GetSpectrum();
            foreach (int scan in spectraData.Keys)
            {
                MS2Spectrum spectrum = spectraData[scan];
                double mz = spectrum.PrecursorMZ();
                int charge = spectrum.PrecursorCharge();

                List<string> candidates = precursorMatch.Match(mz, charge);
                if (candidates.Count == 0)
                    continue;

                List<SearchResult> searched = glycanSearch.Search(spectrum.GetPeaks(), charge, candidates);
                List<SearchResult> results = analyzer.Analyze(searched, mz, scan, 0);
                final.AddRange(results);
                
            }

            //write out
            string outputPath = @"C:\Users\iruiz\Downloads\MSMS\standard.csv";
            //MultiGlycanClassLibrary.util.mass.Glycan.To.SetPermethylation(true, true);
            using (FileStream ostrm = new FileStream(outputPath, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan,glycan,mz,score");
                    foreach (SearchResult r in final.OrderBy(p => p.Scan()))
                    {
                        string output = r.Scan().ToString() + ","
                            + r.Glycan() + ","
                            + r.MZ().ToString() + ","
                            + r.Score().ToString();
                        writer.WriteLine(output);
                    }
                    writer.Flush();
                }
            }

        }

    }
}
