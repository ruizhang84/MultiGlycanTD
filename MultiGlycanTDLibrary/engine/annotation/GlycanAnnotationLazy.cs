using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.engine.search;
using MultiGlycanTDLibrary.model;
using MultiGlycanTDLibrary.model.glycan;
using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;

namespace MultiGlycanTDLibrary.engine.annotation
{
    public class GlycanAnnotationLazy
    {

        GlycanAnnotationSearcher searcher_;
        protected readonly int maxCharge = 3;

        public GlycanAnnotationLazy(GlycanAnnotationSearcher searcher)
        {
            searcher_ = searcher;
        }

        public List<PeakAnnotated> Annotated(List<IPeak> peaks, SearchResult result)
        {
            // init
            searcher_.InitAnnotation(result);

            List<PeakAnnotated> peakAnnotateds = new List<PeakAnnotated>();
            // search peaks
            for (int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                PeakAnnotated annotated = new PeakAnnotated();
                annotated.Peak = peak;
                annotated.Glycan = result.Glycan;

                for (int charge = 1; charge <= Math.Min(maxCharge, result.Charge); charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(), result.Ion, charge);

                    List<GlycanAnnotated> glycans = searcher_.Search(mass);
                    foreach (GlycanAnnotated glycan in glycans)
                    {
                        annotated.Fragments.Add(glycan);
                    }
                }
                if (annotated.Fragments.Count > 0)
                    peakAnnotateds.Add(annotated);
            }

            return peakAnnotateds;
        }
    }
}
