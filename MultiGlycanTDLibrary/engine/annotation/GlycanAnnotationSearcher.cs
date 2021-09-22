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
    public class GlycanAnnotationSearcher
    {
        ISearch<GlycanAnnotated> searcher_;
        Dictionary<string, IGlycan> glycanMaps;
        Dictionary<string, List<Point<GlycanAnnotated>>> annotatedMaps;
        List<FragmentType> types;

        public GlycanAnnotationSearcher(
            ISearch<GlycanAnnotated> searcher,
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
            annotatedMaps = new Dictionary<string, List<Point<GlycanAnnotated>>>();
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
            if (parameter.Permethylated)
            {
                GlycanIonsBuilder.Build.Permethylated = true;
                Glycan.To.SetPermethylation(true, parameter.Reduced);
            }
            else
            {
                GlycanIonsBuilder.Build.Permethylated = false;
                Glycan.To.SetPermethylation(false, parameter.Reduced);
                switch (parameter.Derivatization)
                {
                    case Derivatization.k2AA:
                        GlycanIonsBuilder.Build.Derivatization = GlycanIonsBuilder.k2AA;
                        Glycan.To.Derivatization = Glycan.k2AA;
                        break;
                    case Derivatization.k2AB:
                        GlycanIonsBuilder.Build.Derivatization = GlycanIonsBuilder.k2AB;
                        Glycan.To.Derivatization = Glycan.k2AB;
                        break;
                    default:
                        GlycanIonsBuilder.Build.Derivatization = GlycanIonsBuilder.kWater;
                        Glycan.To.Derivatization = Glycan.kWater;
                        break;
                }
            }





        }

        public void InitAnnotation(SearchResult result)
        {
            string glycan = result.Glycan;
            List<Point<GlycanAnnotated>> points = new List<Point<GlycanAnnotated>>();
            if (annotatedMaps.ContainsKey(glycan))
            {
                points = annotatedMaps[glycan];
            }
            else
            {
                List<double> temp = new();
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
                        temp.Add(mass);
                    }
                }
                annotatedMaps[glycan] = points;
            }
            searcher_.Init(points);
        }

        public List<GlycanAnnotated> Search(double mass)
        {
            return searcher_.SearchContent(mass);
        }

    }
}
