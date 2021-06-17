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
            if (Types.Contains(FragmentTypes.B) || Types.Contains(FragmentTypes.C))
            {
                List<IGlycan> bionsLikeFragments = GlycanFragmentBuilder.BionsLikeFragments(glycan);
                if (Types.Contains(FragmentTypes.B))
                    fragments.AddRange(Bions(bionsLikeFragments));
                if (Types.Contains(FragmentTypes.C))
                    fragments.AddRange(Cions(bionsLikeFragments));
            }

            if (Types.Contains(FragmentTypes.Y) || Types.Contains(FragmentTypes.Z))
            {
                List<IGlycan> yionsLikeFragments = GlycanFragmentBuilder.YionsLikeFragments(glycan);
                if (Types.Contains(FragmentTypes.Y))
                    fragments.AddRange(Yions(yionsLikeFragments));
                if (Types.Contains(FragmentTypes.Z))
                    fragments.AddRange(Zions(yionsLikeFragments));
            }

            if (Types.Contains(FragmentTypes.BY) || Types.Contains(FragmentTypes.BZ) || Types.Contains(FragmentTypes.CY))
            {
                List<IGlycan> byionsLikeFragments = GlycanFragmentBuilder.BYionsLikeFragments(glycan);
                if (Types.Contains(FragmentTypes.BY))
                    fragments.AddRange(BYions(byionsLikeFragments));
                if (Types.Contains(FragmentTypes.BZ))
                    fragments.AddRange(BZions(byionsLikeFragments));
                if (Types.Contains(FragmentTypes.CY))
                    fragments.AddRange(CYions(byionsLikeFragments));
            }

            if (Types.Contains(FragmentTypes.YY) || Types.Contains(FragmentTypes.YZ) || Types.Contains(FragmentTypes.ZZ))
            {
                List<IGlycan> yyionsLikeFragments = GlycanFragmentBuilder.YYionsLikeFragments(glycan);
                if (Types.Contains(FragmentTypes.YY))
                    fragments.AddRange(YYions(yyionsLikeFragments));
                if (Types.Contains(FragmentTypes.YZ))
                    fragments.AddRange(YZions(yyionsLikeFragments));
                if (Types.Contains(FragmentTypes.ZZ))
                    fragments.AddRange(ZZions(yyionsLikeFragments));
            }

            return fragments.Distinct().ToList();
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
