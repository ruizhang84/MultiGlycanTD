using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public class GlycanIonsBuilder
    {
        protected static readonly Lazy<GlycanIonsBuilder>
            lazy = new Lazy<GlycanIonsBuilder>(() => new GlycanIonsBuilder());

        public static GlycanIonsBuilder Build { get { return lazy.Value; } }

        public bool Reduced { get; set; } = true;

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
            List<IGlycan> yionsLikeFragments = GlycanFragmentBuilder.Build.YionsLikeFragments(glycan);
            List<IGlycan> bionsLikeFragments = GlycanFragmentBuilder.Build.BionsLikeFragments(glycan);
            List<IGlycan> yyionsLikeFragments = GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan);
            List<IGlycan> byionsLikeFragments = GlycanFragmentBuilder.Build.BYionsLikeFragments(glycan);

            fragments.AddRange(Bions(bionsLikeFragments));
            fragments.AddRange(Cions(bionsLikeFragments));
            fragments.AddRange(Yions(yionsLikeFragments));
            fragments.AddRange(Zions(yionsLikeFragments));
            fragments.AddRange(BYions(byionsLikeFragments));
            fragments.AddRange(BZions(byionsLikeFragments));
            fragments.AddRange(CYions(byionsLikeFragments));
            fragments.AddRange(YYions(yyionsLikeFragments));
            fragments.AddRange(YZions(yyionsLikeFragments));
            fragments.AddRange(ZZions(yyionsLikeFragments));
            return fragments.Distinct().ToList();
        }

        public List<double> Yions(List<IGlycan> glycans)
        {
            if (Reduced)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced + kHydrogen).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced + kHydrogen).Distinct().ToList(); 
        }

        public List<double> Zions(List<IGlycan> glycans)
        {
            if (Reduced)
                return glycans
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kHydroxyl).Distinct().ToList();
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced - kHydroxyl).Distinct().ToList();
        }

        public List<double> Bions(List<IGlycan> glycans)
        {
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
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) - kWater).Distinct().ToList();
        }

        public List<double> CYions(List<IGlycan> glycans)
        {
            return glycans
                .Select(m => Glycan.To.ComputeFragment(m) + kWater).Distinct().ToList();
        }

        public List<double> CZions(List<IGlycan> glycans)
        {
            return BYions(glycans);
        }
    }
}
