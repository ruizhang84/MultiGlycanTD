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

        ISearch<GlycanAnnotated> searcher_;
        Dictionary<string, IGlycan> glycanMaps;
        List<FragmentType> types;
        protected readonly int maxCharge = 3;

        public GlycanAnnotationLazy(ISearch<GlycanAnnotated> searcher,
            ParameterJson parameter)
        {
            searcher_ = searcher;
            types = new List<FragmentType>()
            {
                FragmentType.B, FragmentType.C, FragmentType.Y, FragmentType.Z,
                FragmentType.BY, FragmentType.BZ, FragmentType.CY,
                FragmentType.YY, FragmentType.YZ, FragmentType.ZZ,
                FragmentType.BYY, FragmentType.BYZ, FragmentType.BZZ,
                FragmentType.CYY, FragmentType.CYZ, FragmentType.CZZ,
                FragmentType.YYY, FragmentType.ZZZ, FragmentType.YYZ, FragmentType.YZZ
            };
            BuilGlycans(parameter);
        }

        public void BuilGlycans(ParameterJson parameter)
        {
            GlycanBuilderSimple glycanBuilder = new GlycanBuilderSimple(
                parameter.HexNAc, parameter.Hex, parameter.Fuc, parameter.NeuAc,
                parameter.NeuGc, parameter.ComplexInclude, parameter.HybridInclude,
                parameter.HighMannoseInclude);
            glycanBuilder.Build();
            glycanMaps = glycanBuilder.GlycanMaps();
        }

        public void InitAnnotation(SearchResult result)
        {

            List<Point<GlycanAnnotated>> points = new List<Point<GlycanAnnotated>>();
            string glycan = result.Glycan;
            foreach (FragmentType type in types)
            {
                List<IGlycan> ionsLikeFragments = GlycanIonsBuilder.Build
                    .FragmentsBuild(glycanMaps[glycan], type);
                foreach (IGlycan g in ionsLikeFragments)
                {
                    double mass = GlycanIonsBuilder.Build.ComputeIonMass(g, type);
                    GlycanAnnotated annoated = new GlycanAnnotated
                    {
                        Parent = glycan,
                        Type = type,
                        Glycan = g.ID()
                    };
                    points.Add(new Point<GlycanAnnotated>(mass, annoated));
                }
            }
            searcher_.Init(points);
        }

        public List<PeakAnnotated> Annotated(List<IPeak> peaks, SearchResult result)
        {
            // init
            InitAnnotation(result);

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

                    List<GlycanAnnotated> glycans = searcher_.SearchContent(mass);
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
