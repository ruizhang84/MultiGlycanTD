using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model.glycan
{
    public class NHighMannose : BaseNGlycan, IGlycan
    {
        //GlcNAc(2) - Man(1) - Fuc - 0, 1, 2,
        // [Man(antenna1) - Man(antenna2)] - 3, 4,
        // [Man(branch1) - Man(branch2) - Man(branch3)]  5 6 7

        private int[] table_ = new int[8];
        
        public override IGlycan Clone()
        {
            IGlycan glycan = new NHighMannose();
            glycan.SetTable(table_);
            glycan.SetComposition(composite);
            return glycan;
        }
        public override bool IsValid()
        {
            // at least three chains
            if (table_[5] == 0 || table_[6] == 0 || table_[7] == 0)
                return false;
            // maksure sorted
            if (table_[5] < table_[6])
                return false;
            return base.IsValid();
        }
        public override int[] Table() { return table_; }
        public override void SetTable(int[] table)
        {
            for (int i = 0; i < table_.Length; i++)
                table_[i] = table[i];
        }

        void AddMonosaccharide(Monosaccharide sugar)
        {
            if (composite.ContainsKey(sugar))
            {
                composite[sugar] += 1;
            }
            else
            {
                composite[sugar] = 1;
            }
        }

        public override List<IGlycan> Grow(Monosaccharide monosaccharide)
        {
            List<IGlycan> glycans = new List<IGlycan>();
            switch (monosaccharide)
            {
                case Monosaccharide.GlcNAc:
                    if (ValidAddGlcNAc())
                    {
                        NHighMannose ptr = CreateByAddGlcNAc();
                        glycans.Add(ptr);
                    }

                    break;
                case Monosaccharide.Man:
                    if (ValidAddManCore())
                    {
                        List<NHighMannose> gs = CreateByAddManCore();
                        glycans.AddRange(gs);
                    }
                    if (ValidAddManBranch())
                    {
                        List<NHighMannose> gs = CreateByAddManBranch();
                        glycans.AddRange(gs);

                    }
                    break;

                case Monosaccharide.Fuc:
                    if (ValidAddFucCore())
                    {
                        NHighMannose ptr = CreateByAddFucCore();
                        glycans.Add(ptr);
                    }
                    break;

                default:
                    break;
            }
            return glycans;
        }


        bool ValidAddGlcNAc()
        {
            return table_[0] < 2;
        }

        NHighMannose CreateByAddGlcNAc()
        {
            var g = new NHighMannose();
            g.SetTable(table_);
            g.table_[0] = g.table_[0] + 1;
            g.SetComposition(composite);
            g.AddMonosaccharide(Monosaccharide.GlcNAc);

            return g;
        }


        bool ValidAddManCore()
        {
            return table_[0] == 2 && (table_[3] == 0 || table_[4] == 0);
        }

        List<NHighMannose> CreateByAddManCore()
        {
            List<NHighMannose> glycans = new List<NHighMannose>();
            if (table_[1] == 0)
            {
                var g = new NHighMannose();
                g.SetTable(table_);
                g.SetComposition(composite);
                g.table_[1] = 1;
                g.AddMonosaccharide(Monosaccharide.Man);
                glycans.Add(g);
            }
            else
            {
                for (int i = 0; i < 2; i++)
                {
                    if (table_[3 + i] == 0)
                    {
                        var g = new NHighMannose();
                        g.SetTable(table_);
                        g.SetComposition(composite);
                        g.table_[3 + i] = 1;
                        g.AddMonosaccharide(Monosaccharide.Man);
                        glycans.Add(g);
                    }
                }
            }
            return glycans;
        }

        bool ValidAddManBranch()
        {
            return table_[3] > 0 || table_[4] > 0;
        }

        List<NHighMannose> CreateByAddManBranch()
        {
            List<NHighMannose> glycans = new List<NHighMannose>();
            for (int i = 0; i < 3; i++)
            {
                if (i == 2 && table_[4] == 0)
                    continue;
                else if (i < 2 && table_[3] == 0)
                    continue;
                var g = new NHighMannose();
                g.SetTable(table_);
                g.table_[i + 5] = g.table_[i + 5] + 1;
                g.SetComposition(composite);
                g.AddMonosaccharide(Monosaccharide.Man);
                glycans.Add(g);  
            }
            return glycans;
        }

        bool ValidAddFucCore()
        {
            return table_[0] >= 1 && table_[2] == 0;  //core
        }

        NHighMannose CreateByAddFucCore()
        {
            var g = new NHighMannose();
            g.SetTable(table_);
            g.table_[2] = 1;
            g.SetComposition(composite);
            g.AddMonosaccharide(Monosaccharide.Fuc);
            return g;
        }

        public override GlycanType Type()
        {
            return GlycanType.NHighMannose;
        }
    }
}
