using MultiGlycanClassLibrary.util.mass;
using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public enum FragmentType
    {
        B, C, Y, Z, BY, BZ, CY, YY, YZ, ZZ,
        BYY, BYZ, BZZ, 
        CYY, CYZ, CZZ,
        YYY, ZZZ, YYZ, YZZ,
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

        public List<FragmentType> Types { get; set; }
            = new List<FragmentType>()
            {
                FragmentType.B, FragmentType.C, FragmentType.Y, FragmentType.Z,
                FragmentType.BY, FragmentType.BZ, FragmentType.CY, FragmentType.YY,
                FragmentType.YZ, FragmentType.ZZ
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
            foreach(FragmentType type in Types)
            {
                fragments.AddRange(Fragments(glycan, type));
            }
            return fragments.Distinct().ToList();
        }

        public List<double> Fragments(IGlycan glycan, FragmentType type)
        {

            if (type == FragmentType.B || type == FragmentType.C)
            {
                List<IGlycan> bionsLikeFragments = GlycanFragmentBuilder.BionsLikeFragments(glycan);
                if (type ==  FragmentType.B)
                    return Bions(bionsLikeFragments);
                if (type == FragmentType.C)
                    return Cions(bionsLikeFragments);
            }

            if (type == FragmentType.Y || type == FragmentType.Z)
            {
                List<IGlycan> yionsLikeFragments = GlycanFragmentBuilder.YionsLikeFragments(glycan);
                if (type == FragmentType.Y)
                    return Yions(yionsLikeFragments);
                if (type == FragmentType.Z)
                    return Zions(yionsLikeFragments);
            }

            if (type == FragmentType.BY || type == FragmentType.BZ || type == FragmentType.CY)
            {
                List<IGlycan> byionsLikeFragments = GlycanFragmentBuilder.BYionsLikeFragments(glycan);
                if (type == FragmentType.BY)
                    return BYions(byionsLikeFragments);
                if (type == FragmentType.BZ)
                    return BZions(byionsLikeFragments);
                if (type == FragmentType.CY)
                    return CYions(byionsLikeFragments);
            }

            if (type == FragmentType.YY || type == FragmentType.YZ || type == FragmentType.ZZ)
            {
                List<IGlycan> yyionsLikeFragments = GlycanFragmentBuilder.YYionsLikeFragments(glycan);
                if (type == FragmentType.YY)
                    return YYions(yyionsLikeFragments);
                if (type == FragmentType.YZ)
                    return YZions(yyionsLikeFragments);
                if (type == FragmentType.ZZ)
                    return ZZions(yyionsLikeFragments);
            }

            if (type == FragmentType.BYY || type == FragmentType.BYZ || type == FragmentType.BZZ ||
                type == FragmentType.CYY || type == FragmentType.CYZ || type == FragmentType.CZZ)
            {
                List<IGlycan> byyionsLikeFragments = GlycanFragmentBuilder.BYYionsLikeFragments(glycan);
                if (type == FragmentType.BYY)
                    return BYYions(byyionsLikeFragments);
                if (type == FragmentType.BYZ)
                    return BYZions(byyionsLikeFragments);
                if (type == FragmentType.BZZ)
                    return BZZions(byyionsLikeFragments);
                if (type == FragmentType.CYY)
                    return BYYions(byyionsLikeFragments);
                if (type == FragmentType.CYZ)
                    return BYZions(byyionsLikeFragments);
                if (type == FragmentType.CZZ)
                    return BZZions(byyionsLikeFragments);
            }

            if (type == FragmentType.YYY || type == FragmentType.YYZ || type == FragmentType.YZZ || 
                type == FragmentType.ZZZ)
            {
                List<IGlycan> yyyionsLikeFragments = GlycanFragmentBuilder.YYYionsLikeFragments(glycan);
                if (type == FragmentType.YYY)
                    return YYYions(yyyionsLikeFragments);
                if (type == FragmentType.YYZ)
                    return YYZions(yyyionsLikeFragments);
                if (type == FragmentType.YZZ)
                    return YZZions(yyyionsLikeFragments);
                if (type == FragmentType.ZZZ)
                    return ZZZions(yyyionsLikeFragments);
            }

            return new List<double>();
        }

        public List<IGlycan> FragmentsBuild(IGlycan glycan, FragmentType type)
        {
            switch (type)
            {
                case FragmentType.B:
                case FragmentType.C:
                    return GlycanFragmentBuilder.BionsLikeFragments(glycan);
                case FragmentType.Y:
                case FragmentType.Z:
                    return GlycanFragmentBuilder.YionsLikeFragments(glycan);
                case FragmentType.BY:
                case FragmentType.BZ:
                case FragmentType.CY:
                    return GlycanFragmentBuilder.BYionsLikeFragments(glycan);
                case FragmentType.YY:
                case FragmentType.YZ:
                case FragmentType.ZZ:
                    return GlycanFragmentBuilder.YYionsLikeFragments(glycan);
                case FragmentType.BYY:
                case FragmentType.BYZ:
                case FragmentType.BZZ:
                case FragmentType.CYY:
                case FragmentType.CYZ:
                case FragmentType.CZZ:
                    return GlycanFragmentBuilder.BYYionsLikeFragments(glycan);
                case FragmentType.YYY:
                case FragmentType.YYZ:
                case FragmentType.YZZ:
                case FragmentType.ZZZ:
                    return GlycanFragmentBuilder.YYYionsLikeFragments(glycan);
            }
            return new List<IGlycan>();
        }

        public double ComputeIonMass(IGlycan subGlycan, FragmentType type)
        {
            switch (type)
            {
                case FragmentType.B:
                    return GlycanIonsBuilder.Build.Bion(subGlycan);
                case FragmentType.C:
                    return GlycanIonsBuilder.Build.Cion(subGlycan);
                case FragmentType.Y:
                    return GlycanIonsBuilder.Build.Yion(subGlycan);
                case FragmentType.Z:
                    return GlycanIonsBuilder.Build.Zion(subGlycan);
                case FragmentType.BY:
                    return GlycanIonsBuilder.Build.BYion(subGlycan);
                case FragmentType.BZ:
                    return GlycanIonsBuilder.Build.BZion(subGlycan);
                case FragmentType.CY:
                    return GlycanIonsBuilder.Build.CYion(subGlycan);
                case FragmentType.YY:
                    return GlycanIonsBuilder.Build.YYion(subGlycan);
                case FragmentType.YZ:
                    return GlycanIonsBuilder.Build.YZion(subGlycan);
                case FragmentType.ZZ:
                    return GlycanIonsBuilder.Build.ZZion(subGlycan);
                case FragmentType.BYY:
                    return GlycanIonsBuilder.Build.BYYion(subGlycan);
                case FragmentType.BYZ:
                    return GlycanIonsBuilder.Build.BYZion(subGlycan);
                case FragmentType.BZZ:
                    return GlycanIonsBuilder.Build.BZZion(subGlycan);
                case FragmentType.CYY:
                    return GlycanIonsBuilder.Build.CYYion(subGlycan);
                case FragmentType.CYZ:
                    return GlycanIonsBuilder.Build.CYZion(subGlycan);
                case FragmentType.CZZ:
                    return GlycanIonsBuilder.Build.CZZion(subGlycan);
                case FragmentType.YYY:
                    return GlycanIonsBuilder.Build.YYYion(subGlycan);
                case FragmentType.YYZ:
                    return GlycanIonsBuilder.Build.YYZion(subGlycan);
                case FragmentType.YZZ:
                    return GlycanIonsBuilder.Build.YZZion(subGlycan);
                case FragmentType.ZZZ:
                    return GlycanIonsBuilder.Build.ZZZion(subGlycan);
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

        // three cleavages
        public double BYYion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan) - kCarbon - kHydrogen * 2;
        }
        public List<double> BYYions(List<IGlycan> glycans)
        {
            return glycans.Select(m => BYYion(m)).Distinct().ToList();
        }

        public double BYZion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan) - kCarbon - kHydrogen * 2 - kWater;
        }

        public List<double> BYZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => BYZion(m)).Distinct().ToList();
        }

        public double BZZion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan) - kWater * 2 - kCarbon - kHydrogen * 2;
        }

        public List<double> BZZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => BZZion(m)).Distinct().ToList();
        }

        public double CYYion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan);
        }

        public List<double> CYYions(List<IGlycan> glycans)
        {
            return glycans.Select(m => CYYion(m)).Distinct().ToList();
        }

        public double CYZion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan);
        }

        public List<double> CYZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => CYZion(m)).Distinct().ToList();
        }

        public double CZZion(IGlycan glycan)
        {
            return Glycan.To.ComputeFragment(glycan);
        }

        public List<double> CZZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => CZZion(m)).Distinct().ToList();
        }

        public double YYYion(IGlycan glycan)
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

        public List<double> YYYions(List<IGlycan> glycans)
        {
            return glycans.Select(m => YYYion(m)).Distinct().ToList();
        }

        public double ZZZion(IGlycan glycan)
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

        public List<double> ZZZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => ZZZion(m)).Distinct().ToList();
        }
        public double YYZion(IGlycan glycan)
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
        public List<double> YYZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => YYZion(m)).Distinct().ToList();
        }

        public double YZZion(IGlycan glycan)
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
        public List<double> YZZions(List<IGlycan> glycans)
        {
            return glycans.Select(m => YZZion(m)).Distinct().ToList();
        }
    }
}
