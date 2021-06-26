﻿using MultiGlycanClassLibrary.util.mass;
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

        // ABEE 149.0841, 2-AA labeled 139.06333, 2-AB labeled 	138.07931
        public double Derivatization { get; set; } = kWater;
        public const double kABEE = 149.0841;
        public const double k2AA = 139.06333;
        public const double k2AB = 138.07931;

        public const double kCarbon = 12.0;
        public const double kOxygen = 15.99491463;
        public const double kHydrogen = 1.007825;
        // methyl
        public const double kMethyl = kCarbon + kHydrogen * 3;
        // oH
        public const double kHydroxyl = kHydrogen + kOxygen;
        // water
        public const double kWater = kHydrogen * 2 + kOxygen;
        // 31
        public const double kNonReduced = kMethyl + kOxygen;
        // 47
        public const double kReduced = kMethyl * 2 + kHydrogen + kOxygen;

        // assume only consist of the part of glycan
        public List<double> Fragments(IGlycan glycan)
        {
            List<double> fragments = new List<double>();
            foreach(FragmentTypes type in Types)
            {
                fragments.AddRange(Fragments(glycan, type));
            }
            return fragments.Distinct().ToList();
        }

        public List<double> Fragments(IGlycan glycan, FragmentTypes type)
        {

            if (type == FragmentTypes.B || type == FragmentTypes.C)
            {
                List<IGlycan> bionsLikeFragments = GlycanFragmentBuilder.BionsLikeFragments(glycan);
                if (type ==  FragmentTypes.B)
                    return Bions(bionsLikeFragments);
                if (type == FragmentTypes.C)
                    return Cions(bionsLikeFragments);
            }

            if (type == FragmentTypes.Y || type == FragmentTypes.Z)
            {
                List<IGlycan> yionsLikeFragments = GlycanFragmentBuilder.YionsLikeFragments(glycan);
                if (type == FragmentTypes.Y)
                    return Yions(yionsLikeFragments);
                if (type == FragmentTypes.Z)
                    return Zions(yionsLikeFragments);
            }

            if (type == FragmentTypes.BY || type == FragmentTypes.BZ || type == FragmentTypes.CY)
            {
                List<IGlycan> byionsLikeFragments = GlycanFragmentBuilder.BYionsLikeFragments(glycan);
                if (type == FragmentTypes.BY)
                    return BYions(byionsLikeFragments);
                if (type == FragmentTypes.BZ)
                    return BZions(byionsLikeFragments);
                if (type == FragmentTypes.CY)
                    return CYions(byionsLikeFragments);
            }

            if (type == FragmentTypes.YY || type == FragmentTypes.YZ || type == FragmentTypes.ZZ)
            {
                List<IGlycan> yyionsLikeFragments = GlycanFragmentBuilder.YYionsLikeFragments(glycan);
                if (type == FragmentTypes.YY)
                    return YYions(yyionsLikeFragments);
                if (type == FragmentTypes.YZ)
                    return YZions(yyionsLikeFragments);
                if (type == FragmentTypes.ZZ)
                    return ZZions(yyionsLikeFragments);
            }

            return new List<double>();
        }

        public List<IGlycan> FragmentsBuild(IGlycan glycan, FragmentTypes type)
        {
            switch (type)
            {
                case FragmentTypes.B:
                case FragmentTypes.C:
                    return GlycanFragmentBuilder.BionsLikeFragments(glycan);
                case FragmentTypes.Y:
                case FragmentTypes.Z:
                    return GlycanFragmentBuilder.YionsLikeFragments(glycan);
                case FragmentTypes.BY:
                case FragmentTypes.BZ:
                case FragmentTypes.CY:
                    return GlycanFragmentBuilder.BYionsLikeFragments(glycan);
                case FragmentTypes.YY:
                case FragmentTypes.YZ:
                case FragmentTypes.ZZ:
                    return GlycanFragmentBuilder.YYionsLikeFragments(glycan);

            }
            return new List<IGlycan>();
        }

        public double ComputeIonMass(IGlycan subGlycan, FragmentTypes type)
        {
            switch (type)
            {
                case FragmentTypes.B:
                    return GlycanIonsBuilder.Build.Bion(subGlycan);
                case FragmentTypes.C:
                    return GlycanIonsBuilder.Build.Cion(subGlycan);
                case FragmentTypes.Y:
                    return GlycanIonsBuilder.Build.Yion(subGlycan);
                case FragmentTypes.Z:
                    return GlycanIonsBuilder.Build.Zion(subGlycan);
                case FragmentTypes.BY:
                    return GlycanIonsBuilder.Build.BYion(subGlycan);
                case FragmentTypes.BZ:
                    return GlycanIonsBuilder.Build.BZion(subGlycan);
                case FragmentTypes.CY:
                    return GlycanIonsBuilder.Build.CYion(subGlycan);
                case FragmentTypes.YY:
                    return GlycanIonsBuilder.Build.YYion(subGlycan);
                case FragmentTypes.YZ:
                    return GlycanIonsBuilder.Build.YZion(subGlycan);
                case FragmentTypes.ZZ:
                    return GlycanIonsBuilder.Build.ZZion(subGlycan);
            }
            return 0;
        }

        public double Yion(IGlycan glycan)
        {
            if (Permethylated)
            {
                if(Reduced)
                {
                    return Glycan.To.ComputeFragment(glycan) + kReduced + kHydrogen;
                }
                else
                {
                    return Glycan.To.ComputeFragment(glycan) + kNonReduced + kHydrogen;
                }
            }
            return Glycan.To.ComputeFragment(glycan) + Derivatization;
        }

        public List<double> Yions(List<IGlycan> glycans)
        {
            return glycans.Select(m => Yion(m)).Distinct().ToList();
        }

        public double Zion(IGlycan glycan)
        {
            if (Permethylated)
            {
                if (Reduced)
                {
                    return Glycan.To.ComputeFragment(glycan) + kReduced - kHydroxyl;
                }
                else
                {
                    return Glycan.To.ComputeFragment(glycan) + kNonReduced - kHydroxyl;
                }
            }
            return Glycan.To.ComputeFragment(glycan) + Derivatization - kWater;
        }

        public List<double> Zions(List<IGlycan> glycans)
        {
            return glycans.Select(m => Zion(m)).Distinct().ToList();
        }

        public double Bion(IGlycan glycan)
        {
            if (Permethylated)
            {
                return Glycan.To.ComputeFragment(glycan) + kCarbon + kHydrogen * 2;
            }
            return Glycan.To.ComputeFragment(glycan);
        }


        public List<double> Bions(List<IGlycan> glycans)
        {
            return glycans.Select(m => Bion(m)).Distinct().ToList();
        }

        public double Cion(IGlycan glycan)
        {
            if (Permethylated)
                return Glycan.To.ComputeFragment(glycan) + kCarbon + kHydrogen * 4 + kOxygen;
            return Glycan.To.ComputeFragment(glycan) + kWater;
        }


        public List<double> Cions(List<IGlycan> glycans)
        {
            return glycans.Select(m => Cion(m)).Distinct().ToList();
        }

        public double YYion(IGlycan glycan)
        {
            if (Permethylated)
            {
                if (Reduced)
                {
                    return Glycan.To.ComputeFragment(glycan) + kReduced - kCarbon - kHydrogen;
                }
                else
                {
                    return Glycan.To.ComputeFragment(glycan) + kNonReduced - kCarbon - kHydrogen;
                }
            }
            return Glycan.To.ComputeFragment(glycan) + Derivatization;
        }

        public List<double> YYions(List<IGlycan> glycans)
        {
            return glycans.Select(m => YYion(m)).Distinct().ToList();
        }

        public double YZion(IGlycan glycan)
        {
            if (Permethylated)
            {
                if (Reduced)
                {
                    return Glycan.To.ComputeFragment(glycan) + kReduced - kOxygen - kCarbon - kHydrogen * 3;
                }
                else
                {
                    return Glycan.To.ComputeFragment(glycan);
                }
            }
            return Glycan.To.ComputeFragment(glycan) + Derivatization - kWater;
        }

        public List<double> YZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => YZion(m)).Distinct().ToList(); ;
        }

        public double ZZion(IGlycan glycan)
        {
            if (Permethylated)
            {
                if (Reduced)
                {
                    return Glycan.To.ComputeFragment(glycan) + kReduced - kOxygen * 2 - kCarbon - kHydrogen * 5;
                }
                else
                {
                    return Glycan.To.ComputeFragment(glycan) + kNonReduced - kOxygen * 2 - kCarbon - kHydrogen * 5;
                }
            }
            return Glycan.To.ComputeFragment(glycan) + Derivatization - kWater * 2;            
        }

        public List<double> ZZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => ZZion(m)).Distinct().ToList();
        }

        public double BYion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan);
        }

        public List<double> BYions(List<IGlycan> glycans)
        {
            return glycans.Select(m => BYion(m)).Distinct().ToList(); ;
        }

        public double BZion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan) - kWater;
        }


        public List<double> BZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => BZion(m)).Distinct().ToList(); ;
        }

        public double CYion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan) + kWater;       
        }


        public List<double> CYions(List<IGlycan> glycans)
        {
            return glycans.Select(m => CYion(m)).Distinct().ToList();
        }

        public double CZion(IGlycan glycan)
        {
            return BYion(glycan);
        }

        public List<double> CZions(List<IGlycan> glycans)
        {
            return BYions(glycans);
        }
    }
}
