using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.model;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.annotation
{
    public class GlycanAnnotation
    {
        ISearch<GlycanAnnotated> searcher_;
        //private readonly int maxCharge = 2;
        private Random random;
        private readonly int lower = 1;
        private readonly int upper = 30;
        public int Seed { get; set; } = 2;

        public GlycanAnnotation(ISearch<GlycanAnnotated> searcher, 
            Dictionary<double, List<GlycanAnnotated>> massMap)
        {
            searcher_ = searcher;
            List<Point<GlycanAnnotated>> points = new List<Point<GlycanAnnotated>>();
            foreach (double mass in massMap.Keys)
            {
                foreach (GlycanAnnotated glycan in massMap[mass])
                {
                    points.Add(new Point<GlycanAnnotated>(mass, glycan));
                }
            }
            searcher_.Init(points);
            random = new Random(Seed);
        }

        public List<PeakAnnotated> Annotated(List<IPeak> peaks, int precursorCharge,
            List<SearchResult> candidates, bool decoy = false, double ion = 1.0078)
        {
            // init
            Dictionary<string, PeakAnnotated> results = new Dictionary<string, PeakAnnotated>();
            HashSet<string> glycanCandid = new HashSet<string>();
            foreach (SearchResult r in candidates)
            {
                glycanCandid.UnionWith(r.Isomers());
            }
            // search peaks
            for (int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                double randomMass = 0;
                if (decoy)
                    randomMass = random.NextDouble() * (upper - lower) + lower;
                for (int charge = 1; charge <= precursorCharge; charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(), ion, charge);
                    if (decoy)
                        mass += randomMass;

                    List<GlycanAnnotated> glycans = searcher_.Search(mass);
                    foreach (GlycanAnnotated glycan in glycans)
                    {
                        if (glycanCandid.Contains(glycan.Parent))
                        {
                            if (!results.ContainsKey(glycan.Parent))
                            {
                                results[glycan.Parent] = new PeakAnnotated();
                                results[glycan.Parent].Peak = peak;
                                results[glycan.Parent].Glycan = glycan.Parent;
                            }
                            results[glycan.Parent].Fragments.Add(glycan);
                        }
                    }
                }

            }
           
            return results.Values.ToList();
        }
    }
}
