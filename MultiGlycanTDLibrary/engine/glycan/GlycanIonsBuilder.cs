using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public enum FragmentTypes
    {
        B, C, Y, Z, BY, BZ, CY, YY, YZ, ZZ
    }


    public class GlycanIonsBuilder
    {
        protected static readonly Lazy<GlycanIonsBuilder>
            lazy = new Lazy<GlycanIonsBuilder>(() => new GlycanIonsBuilder());

        public static GlycanIonsBuilder Build { get { return lazy.Value; } }

        public Dictionary<string, List<double>> FragmentsMap { get; set; }
            = new Dictionary<string, List<double>>();

        public bool Reduced { get; set; } = true;
        public bool Permethylated { get; set; } = true;
        public List<FragmentTypes> Types { get; set; }
            = new List<FragmentTypes>()
            {
                FragmentTypes.B, FragmentTypes.C, FragmentTypes.Y, FragmentTypes.Z,
                FragmentTypes.BY, FragmentTypes.BZ, FragmentTypes.CY, FragmentTypes.YY,
                FragmentTypes.YZ, FragmentTypes.ZZ
            };


        object obj = new object();

        public const double kCarbon = 12.0;
        public const double kOxygen = 15.99491463;
        public const double kHydrogen = 1.007825;
        // methyl
        public const double kMethyl = kCarbon + kHydrogen * 3;
        // H2O
        public const double kWater = kHydrogen * 2 + kOxygen;
        // oH
        public const double kHydroxyl = kHydrogen + kOxygen;
        // 31
        public const double kNonReduced = kMethyl + kOxygen;
        // 47
        public const double kReduced = kMethyl * 2 + kHydrogen + kOxygen;

        // assume only consist of the part of glycan
        public List<double> Fragments(IGlycan glycan)
        {
            List<double> fragments = new List<double>();
            if (Types.Contains(FragmentTypes.B) || Types.Contains(FragmentTypes.C))
            {
                List<IGlycan> bionsLikeFragments = GlycanFragmentBuilder.Build.BionsLikeFragments(glycan);
                if (Types.Contains(FragmentTypes.B))
                    fragments.AddRange(Bions(bionsLikeFragments));
                if (Types.Contains(FragmentTypes.C))
                    fragments.AddRange(Cions(bionsLikeFragments));
            }

            if (Types.Contains(FragmentTypes.Y) || Types.Contains(FragmentTypes.Z))
            {
                List<IGlycan> yionsLikeFragments = GlycanFragmentBuilder.Build.YionsLikeFragments(glycan);
                if (Types.Contains(FragmentTypes.Y))
                    fragments.AddRange(Yions(yionsLikeFragments));
                if (Types.Contains(FragmentTypes.Z))
                    fragments.AddRange(Zions(yionsLikeFragments));
            }

            if (Types.Contains(FragmentTypes.BY) || Types.Contains(FragmentTypes.BZ) || Types.Contains(FragmentTypes.CY))
            {
                List<IGlycan> byionsLikeFragments = GlycanFragmentBuilder.Build.BYionsLikeFragments(glycan);
                if (Types.Contains(FragmentTypes.BY))
                    fragments.AddRange(BYions(byionsLikeFragments));
                if (Types.Contains(FragmentTypes.BZ))
                    fragments.AddRange(BZions(byionsLikeFragments));
                if (Types.Contains(FragmentTypes.CY))
                    fragments.AddRange(CYions(byionsLikeFragments));
            }

            if (Types.Contains(FragmentTypes.YY) || Types.Contains(FragmentTypes.YZ) || Types.Contains(FragmentTypes.ZZ))
            {
                List<IGlycan> yyionsLikeFragments = GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan);
                if (Types.Contains(FragmentTypes.YY))
                    fragments.AddRange(YYions(yyionsLikeFragments));
                if (Types.Contains(FragmentTypes.YZ))
                    fragments.AddRange(YZions(yyionsLikeFragments));
                if (Types.Contains(FragmentTypes.ZZ))
                    fragments.AddRange(ZZions(yyionsLikeFragments));
            }

            return fragments.Distinct().ToList();
        }

        public List<double> Yions(List<IGlycan> glycans)
        {
            if (!Permethylated)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
            if (Reduced)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced + kHydrogen).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced + kHydrogen).Distinct().ToList(); 
        }

        public List<double> Zions(List<IGlycan> glycans)
        {
            if (!Permethylated)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
            if (Reduced)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kHydroxyl).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced - kHydroxyl).Distinct().ToList();
        }

        public List<double> Bions(List<IGlycan> glycans)
        {
            if (!Permethylated)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kCarbon + kHydrogen * 2).Distinct().ToList();
        }

        public List<double> Cions(List<IGlycan> glycans)
        {
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kCarbon + kHydrogen * 4 + kOxygen).Distinct().ToList();
        }

        public List<double> YYions(List<IGlycan> glycans)
        {
            if (!Permethylated)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
            if (Reduced)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kCarbon - kHydrogen).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced - kCarbon - kHydrogen).Distinct().ToList();
        }

        public List<double> YZions(List<IGlycan> glycans)
        {
            if (Reduced)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kOxygen - kCarbon - kHydrogen * 3).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
        }

        public List<double> ZZions(List<IGlycan> glycans)
        {
            if (!Permethylated)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
            if (Reduced)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kOxygen * 2 - kCarbon - kHydrogen * 5).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced - kOxygen * 2 - kCarbon - kHydrogen * 5).Distinct().ToList();
        }

        public List<double> BYions(List<IGlycan> glycans)
        {
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
        }

        public List<double> BZions(List<IGlycan> glycans)
        {
            if (!Permethylated)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) - kWater).Distinct().ToList();
        }

        public List<double> CYions(List<IGlycan> glycans)
        {
            if (!Permethylated)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kWater).Distinct().ToList();
        }

        public List<double> CZions(List<IGlycan> glycans)
        {
            return BYions(glycans);
        }
    }
}
