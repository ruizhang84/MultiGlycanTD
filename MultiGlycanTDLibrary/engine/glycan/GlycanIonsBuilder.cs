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
            return GlycanFragmentBuilder.Build.YionsLikeFragments(glycan)
                .Select(m => m.Mass() + kHydrogen).Distinct().ToList();
        }

        public List<double> Zions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.YionsLikeFragments(glycan)
                .Select(m => m.Mass() - kHydroxyl).Distinct().ToList();
        }

        public List<double> Bions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BionsLikeFragments(glycan)
                .Select(m => m.Mass() + kCarbon + kHydrogen * 2).Distinct().ToList();
        }

        public List<double> Cions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BionsLikeFragments(glycan)
                .Select(m => m.Mass() + kCarbon + kHydrogen * 4 + kOxygen).Distinct().ToList();
        }

        public List<double> YYions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                .Select(m => m.Mass() - kCarbon - kHydrogen).Distinct().ToList();
        }

        public List<double> YZions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                .Select(m => m.Mass() - kOxygen - kCarbon - kHydrogen * 3).Distinct().ToList();
        }

        public List<double> ZZions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.YYionsLikeFragments(glycan)
                .Select(m => m.Mass() - kReduced - kHydrogen * 2).Distinct().ToList();
        }

        public List<double> BYions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BYionsLikeFragments(glycan)
                .Select(m => m.Mass()).Distinct().ToList();
        }

        public List<double> BZions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BYionsLikeFragments(glycan)
                .Select(m => m.Mass() - kWater).Distinct().ToList();
        }

        public List<double> CYions(IGlycan glycan)
        {
            return GlycanFragmentBuilder.Build.BYionsLikeFragments(glycan)
                .Select(m => m.Mass() + kWater).Distinct().ToList();
        }

        public List<double> CZions(IGlycan glycan)
        {
            return BYions(glycan);
        }


       


    }
}
