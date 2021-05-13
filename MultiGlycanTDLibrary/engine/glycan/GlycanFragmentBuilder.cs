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

        public static GlycanFragmentBuilder Build { get { return lazy.Value; } }

        private int CountYFragments(IGlycan sub, int[] table, int limit)
        {
            int diff = 0;
            int[] subTable = sub.Table();
            // make sure not fragments on core, GlcNAc(2) - Man(3)
            if (subTable[0] + subTable[1] < 5)
                return 0;

            switch (sub.Type())
            {
                case GlycanType.NGlycanComplex:
                    for (int i = 0; i < 4; i++)
                    {
                        // branch
                        if (table[i + 4] != subTable[i + 4] ||
                            table[i + 8] != subTable[i + 8] ||
                            table[i + 16] != subTable[i + 16] ||
                            table[i + 20] != subTable[i + 20])
                            diff++;
                        // fuc count as additional
                        if (table[i + 12] != subTable[i + 12])
                            diff++;
                        if (diff > limit)
                            break;
                    }
                    break;
                case GlycanType.NGlycanHybrid:
                    for (int i = 0; i < 2; i++)
                    {
                        if (table[i + 6] != subTable[i + 6] ||
                            table[i + 8] != subTable[i + 8] ||
                            table[i + 12] != subTable[i + 12] ||
                            table[i + 14] != subTable[i + 14])
                            diff++;
                        if (table[i + 4] != subTable[i + 4])
                            diff++;
                        if (table[i + 10] != subTable[i + 11])
                            diff++;
                        if (diff > limit)
                            break;
                    }
                    break;
                case GlycanType.NHighMannose:
                    for (int i = 0; i < 3; i++)
                    {
                        if (table[i + 3] != subTable[i + 3])
                            diff++;
                        if (diff > limit)
                            break;
                    }
                    break;
            }
            return diff;
        }

        public List<IGlycan> YionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanYFragment = new List<IGlycan>();
            int[] table = glycan.Table();
            foreach (IGlycan sub in glycan.Children())
            {
                int diff = CountYFragments(sub, table, 1);
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
            int[] table = glycan.Table();
            foreach (IGlycan sub in glycan.Children())
            {
                int diff = CountYFragments(sub, table, 2);
                if (diff == 2)
                {
                    glycanYFragment.Add(sub);
                }
            }
            return glycanYFragment;
        }

        private List<IGlycan> BionsBranchFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBFragment = new List<IGlycan>();
            int[] table = glycan.Table();
            IGlycan newGlycan;
            int[] subTable;
            SortedDictionary<Monosaccharide, int> compose;

            // cleave on terminal - single branch
            switch (glycan.Type())
            {
                case GlycanType.NGlycanComplex:
                    // branch exists
                    if (table[4] == 0)
                        break;
                    // compose
                    compose = new SortedDictionary<Monosaccharide, int>();
                    compose[Monosaccharide.GlcNAc] = table[0] + table[3] + table[4];
                    compose[Monosaccharide.Gal] = table[8];
                    compose[Monosaccharide.Man] = table[1];
                    compose[Monosaccharide.Fuc] = table[2] + table[12];
                    compose[Monosaccharide.NeuAc] = table[16];
                    compose[Monosaccharide.NeuGc] = table[20];
                    // core
                    subTable = new int[table.Length];
                    for (int i = 0; i < 4; i++)
                    {
                        subTable[i] = table[i];
                    }
                    // longest branch
                    subTable[4] = table[4];
                    subTable[8] = table[8];
                    subTable[12] = table[12];
                    subTable[16] = table[16];
                    subTable[20] = table[20];
                    // set glycan
                    newGlycan = glycan.Clone();
                    newGlycan.SetTable(subTable);
                    newGlycan.SetComposition(compose);
                    glycanBFragment.Add(newGlycan);
                    break;
                case GlycanType.NGlycanHybrid:
                    if (table[6] == 0)
                        break;
                    // core
                    subTable = new int[table.Length];
                    for (int i = 0; i < 4; i++)
                    {
                        subTable[i] = table[i];
                    }
                    // compose
                    compose = new SortedDictionary<Monosaccharide, int>();
                    compose[Monosaccharide.GlcNAc] = table[0] + table[3] + table[6];
                    compose[Monosaccharide.Gal] = table[8];
                    compose[Monosaccharide.Man] = table[1];
                    compose[Monosaccharide.Fuc] = table[2] + table[10];
                    compose[Monosaccharide.NeuAc] = table[12];
                    compose[Monosaccharide.NeuGc] = table[14];
                    // longest branch
                    subTable[6] = table[6];
                    subTable[8] = table[8];
                    subTable[10] = table[10];
                    subTable[12] = table[12];
                    subTable[14] = table[14];
                    // set
                    newGlycan = glycan.Clone();
                    newGlycan.SetTable(subTable);
                    newGlycan.SetComposition(compose);
                    glycanBFragment.Add(newGlycan);
                    break;
                case GlycanType.NHighMannose:
                    // check brach exists
                    if (table[3] == 0)
                        break;
                    subTable = new int[table.Length];
                    compose = new SortedDictionary<Monosaccharide, int>();
                    compose[Monosaccharide.GlcNAc] = table[0];
                    compose[Monosaccharide.Man] = table[1] + table[3];
                    compose[Monosaccharide.Fuc] = table[2];
                    for (int i = 0; i < 4; i++)
                    {
                        subTable[i] = table[i];
                    }
                    newGlycan = glycan.Clone();
                    newGlycan.SetTable(subTable);
                    newGlycan.SetComposition(compose);
                    glycanBFragment.Add(newGlycan);
                    break;
            }
            return glycanBFragment;
        }

        private List<IGlycan> BionsCoreFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBFragment = new List<IGlycan>();
            // GlcNAc (Fuc) + Man + Man(2)
            IGlycan newGlycan = glycan.Clone();
            SortedDictionary<Monosaccharide, int> compose = newGlycan.Composition();
            int[] subTable = newGlycan.Table();
            compose[Monosaccharide.GlcNAc] -= 1;
            subTable[0] -= 1;
            newGlycan.SetTable(subTable);
            newGlycan.SetComposition(compose);
            glycanBFragment.Add(newGlycan);

            return glycanBFragment;
        }

        public List<IGlycan> BionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBFragment = new List<IGlycan>();
            int[] table = glycan.Table();
            // check branch exists
            if (table[0] + table[1] < 5)
                return glycanBFragment;
            // cleave on core
            glycanBFragment.AddRange(BionsCoreFragments(glycan));

            // cleave on branch
            glycanBFragment.AddRange(BionsBranchFragments(glycan));  
            return glycanBFragment;
        }

        public List<IGlycan> BYionsLikeFragments(IGlycan glycan)
        {
            List<IGlycan> glycanBYFragment = new List<IGlycan>();

            List<IGlycan> YionsFragments = YionsLikeFragments(glycan);
            foreach(IGlycan sub in YionsFragments)
            {
                glycanBYFragment.AddRange(BionsCoreFragments(sub));

            }
            return glycanBYFragment;
        }



    }
}
