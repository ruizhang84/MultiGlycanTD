using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public class GlycanFragmentBuilder
    {
        public static List<IGlycan> YionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanYFragment = new List<IGlycan>();       
            foreach (IGlycan sub in glycan.Fragments())
            {
                int diff = GlycanFragmentBuilderHelper.CountYCut(sub, glycan, 1);
                if (diff == 1)
                {
                    glycanYFragment.Add(sub);
                }
                
            }
            return glycanYFragment;
        }
        public static List<IGlycan> YYionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanYFragment = new List<IGlycan>();
            foreach (IGlycan sub in glycan.Fragments())
            {
                int diff = GlycanFragmentBuilderHelper.CountYCut(sub, glycan, 2);
                if (diff == 2)
                {
                    glycanYFragment.Add(sub);
                }
            }
            return glycanYFragment;
        }
        public static List<IGlycan> YYYionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanYFragment = new List<IGlycan>();
            foreach (IGlycan sub in glycan.Fragments())
            {
                int diff = GlycanFragmentBuilderHelper.CountYCut(sub, glycan, 3);
                if (diff == 3)
                {
                    glycanYFragment.Add(sub);
                }
            }
            return glycanYFragment;
        }
        public static List<IGlycan> BionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBFragment = new List<IGlycan>();
            List<IGlycan> yionsFragments = YionsLikeFragments(glycan);
            foreach (IGlycan sub in yionsFragments)
            {
                glycanBFragment.Add(GlycanFragmentBuilderHelper.ComplementaryFragment(sub, glycan));
            }
            return glycanBFragment;
        }
        public static List<IGlycan> BYionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBYFragment = new List<IGlycan>();

            List<IGlycan> YionsFragments = YionsLikeFragments(glycan);
            foreach(IGlycan sub in YionsFragments)
            {
                List<IGlycan> subYionsFragments = YionsLikeFragments(sub);
                foreach(IGlycan subSub in subYionsFragments)
                {
                    if(GlycanFragmentBuilderHelper.ContainsCut(glycan, sub, subSub))
                    {
                        glycanBYFragment.Add(
                            GlycanFragmentBuilderHelper.ComplementaryFragment(subSub, sub));
                    }
                }
            }
            return glycanBYFragment;
        }
        public static List<IGlycan> BYYionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBYFragment = new List<IGlycan>();

            List<IGlycan> BYionsFragments = BYionsLikeFragments(glycan);
            foreach (IGlycan sub in BYionsFragments)
            {
                List<IGlycan> subYionsFragments = YionsLikeFragments(sub);
                glycanBYFragment.AddRange(subYionsFragments);
            }
            return glycanBYFragment;
        }


    }
}
