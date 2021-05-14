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
        // assume N-Glycan with pentacore
        protected static readonly Lazy<GlycanFragmentBuilder>
              lazy = new Lazy<GlycanFragmentBuilder>(() => new GlycanFragmentBuilder());

        // fragments higher than 3. 
        public static GlycanFragmentBuilder Build { get { return lazy.Value; } }

        public List<IGlycan> YionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanYFragment = new List<IGlycan>();

            foreach (IGlycan sub in glycan.Children())
            {
                int diff = GlycanFragmentBuilderHelper.CountYCut(sub, glycan, 1);
                if (diff == 1)
                {
                    glycanYFragment.Add(sub);
                }
            }
            return glycanYFragment;
        }
        public List<IGlycan> YYionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanYFragment = new List<IGlycan>();

            foreach (IGlycan sub in glycan.Children())
            {
                int diff = GlycanFragmentBuilderHelper.CountYCut(sub, glycan, 2);
                if (diff == 2)
                {
                    glycanYFragment.Add(sub);
                }
            }
            return glycanYFragment;
        }
        public List<IGlycan> BionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBFragment = new List<IGlycan>();
            List<IGlycan> yionsFragments = YionsLikeFragments(glycan);
            foreach (IGlycan sub in yionsFragments)
            {
                glycanBFragment.Add(GlycanFragmentBuilderHelper.ComplementaryFragment(sub, glycan));
            }
            return glycanBFragment;
        }
        public List<IGlycan> BYionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBYFragment = new List<IGlycan>();

            List<IGlycan> YionsFragments = YionsLikeFragments(glycan);
            foreach(IGlycan sub in YionsFragments)
            {
                glycanBYFragment.AddRange(BionsLikeFragments(sub));

            }
            return glycanBYFragment;
        }



    }
}
