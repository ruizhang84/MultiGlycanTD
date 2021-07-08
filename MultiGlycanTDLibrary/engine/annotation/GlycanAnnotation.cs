using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MultiGlycanTDLibrary.engine.annotation
{
    public class GlycanAnnotation
    {
        ISearch<GlycanAnnotated> searcher_;
        private readonly int maxCharge = 3;
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
        }

        public List<PeakAnnotated> Annotated(List<IPeak> peaks, int precursorCharge,
            List<SearchResult> candidates, double ion = 1.0078)
        {
            // init
            List<PeakAnnotated> peakAnnotateds = new List<PeakAnnotated>();
            List<string> glycanCandid = new List<string>();
            foreach (SearchResult r in candidates)
            {
                glycanCandid.Add(r.Glycan);
            }
            // search peaks
            for (int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                Dictionary<string, PeakAnnotated> results = new Dictionary<string, PeakAnnotated>();
                for (int charge = 1; charge <= Math.Min(maxCharge, precursorCharge); charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(), ion, charge);

                    List<GlycanAnnotated> glycans = searcher_.SearchContent(mass);
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
                peakAnnotateds.AddRange(results.Values.ToList());

            }

            return peakAnnotateds;
        }
    }
}
