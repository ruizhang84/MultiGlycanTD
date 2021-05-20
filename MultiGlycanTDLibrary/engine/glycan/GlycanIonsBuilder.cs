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
            fragments.AddRange(Bions(glycan));
            fragments.AddRange(Cions(glycan));
            fragments.AddRange(Yions(glycan));
            fragments.AddRange(Zions(glycan));
            fragments.AddRange(BYions(glycan));
            fragments.AddRange(BZions(glycan));
            fragments.AddRange(CYions(glycan));
            fragments.AddRange(YYions(glycan));
            fragments.AddRange(YZions(glycan));
            fragments.AddRange(ZZions(glycan));
            return fragments.Distinct().ToList();
        }

        public List<double> ExtraFragments(IGlycan glycan)
        {
            List<double> fragments = new List<double>();
            fragments.AddRange(Yions(glycan));
            return fragments.Distinct().ToList();
        }

        public List<double> Yions(IGlycan glycan)
        {
            if (Reduced)
                return GlycanFragmentBuilder.Build.YionsLikeFragments(glycan)
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced + kHydrogen).Distinct().ToList();
            return GlycanFragmentBuilder.Build.YionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced + kHydrogen).Distinct().ToList(); 
        }

        public List<double> Zions(IGlycan glycan)
        {
            if (Reduced)
                return GlycanFragmentBuilder.Build.YionsLikeFragments(glycan)
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kHydroxyl).Distinct().ToList();
            return GlycanFragmentBuilder.Build.YionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced - kHydroxyl).Distinct().ToList();
        }

        public List<double> Bions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m) + kCarbon + kHydrogen * 2).Distinct().ToList();
        }

        public List<double> Cions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m) + kCarbon + kHydrogen * 4 + kOxygen).Distinct().ToList();
        }

        public List<double> YYions(IGlycan glycan)
        {
            if (Reduced)
                return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kCarbon - kHydrogen).Distinct().ToList();
            return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced - kCarbon - kHydrogen).Distinct().ToList();
        }

        public List<double> YZions(IGlycan glycan)
        {
            if (Reduced)
                return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kOxygen - kCarbon - kHydrogen * 3).Distinct().ToList();
            return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
        }

        public List<double> ZZions(IGlycan glycan)
        {
            if (Reduced)
                return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                    .Select(m => Glycan.To.ComputeFragment(m) + kReduced - kOxygen * 2 - kCarbon - kHydrogen * 5).Distinct().ToList();
            return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m) + kNonReduced - kOxygen * 2 - kCarbon - kHydrogen * 5).Distinct().ToList();
        }

        public List<double> BYions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BYionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m)).Distinct().ToList();
        }

        public List<double> BZions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BYionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m) - kWater).Distinct().ToList();
        }

        public List<double> CYions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BYionsLikeFragments(glycan)
                .Select(m => Glycan.To.ComputeFragment(m) + kWater).Distinct().ToList();
        }

        public List<double> CZions(IGlycan glycan)
        {
            return BYions(glycan);
        }
    }
}
