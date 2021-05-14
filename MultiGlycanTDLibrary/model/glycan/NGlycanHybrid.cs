using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model.glycan
{
    public class NGlycanHybrid : BaseNGlycan, IGlycan
    {
        //GlcNAc(2) - Man(3) - Fuc(1) - GlcNAc(bisect,1) 0, 1, 2, 3
        //[Man(antenna1) - Man(antenna2)] - 4, 5, 
        //[Man(branch1) - Man(branch2)]  6, 7
        //[GlcNAc(branch1) - GlcNAc(branch2)] 8, 9
        //[Gal(branch1) - Gal(branch2)]   10, 11
        //[Fuc, Fuc]         12, 13
        //[NeuAc, NeuAc]     14, 15
        //[NeuGc, NeuGc]     16, 17

        private int[] table_ = new int[18];

        public override IGlycan Clone()
        {
            IGlycan glycan = new NGlycanHybrid();
            glycan.SetTable(table_);
            glycan.SetComposition(composite);
            return glycan;
        }
        public override bool IsValid()
        {
            // at least two chains
            if (table_[6] == 0 || table_[8] == 0)
                return false;
            // maksure sorted
            if (table_[4] < table_[5])
                return false;
            if (table_[6] < table_[7])
                return false;
            if (table_[8] < table_[8])
                return false;
            if (table_[10] < table_[11])
                return false;
            if (table_[12] < table_[13])
                return false;
            if (table_[14] < table_[15])
                return false;
            if (table_[16] < table_[17])
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
                    if (ValidAddGlcNAcCore())
                    {
                        NGlycanHybrid ptr = CreateByAddGlcNAcCore();
                        glycans.Add(ptr);
                    }
                    else if (ValidAddGlcNAc())
                    {
                        if (ValidAddGlcNAcBisect())
                        {
                            NGlycanHybrid ptr = CreateByAddGlcNAcBisect();
                            glycans.Add(ptr);
                        }
                        if (ValidAddGlcNAcBranch())
                        {
                            List<NGlycanHybrid> gs = CreateByAddGlcNAcBranch();
                            glycans.AddRange(gs);
                        }
                    }
                    break;

                case Monosaccharide.Man:
                    if (ValidAddManCore())
                    {
                        List<NGlycanHybrid> gs = CreateByAddManCore();
                        glycans.AddRange(gs);
                    }
                    else if (ValidAddManBranch())
                    {
                        List<NGlycanHybrid> gs = CreateByAddManBranch();
                        glycans.AddRange(gs);

                    }
                    break;

                case Monosaccharide.Gal:
                    if (ValidAddGal())
                    {
                        List<NGlycanHybrid> gs = CreateByAddGal();
                        glycans.AddRange(gs);
                    }
                    break;

                case Monosaccharide.Fuc:
                    if (ValidAddFucCore())
                    {
                        NGlycanHybrid ptr = CreateByAddFucCore();
                        glycans.Add(ptr);
                    }
                    else if (ValidAddFucTerminal())
                    {
                        List<NGlycanHybrid> gs = CreateByAddFucTerminal();
                        glycans.AddRange(gs);
                    }
                    break;

                case Monosaccharide.NeuAc:
                    if (ValidAddNeuAc())
                    {
                        List<NGlycanHybrid> gs = CreateByAddNeuAc();
                        glycans.AddRange(gs);
                    }
                    break;

                case Monosaccharide.NeuGc:
                    if (ValidAddNeuGc())
                    {
                        List<NGlycanHybrid> gs = CreateByAddNeuGc();
                        glycans.AddRange(gs);
                    }
                    break;

                default:
                    break;
            }
            return glycans;
        }


        bool ValidAddGlcNAcCore()
        {
            return table_[0] < 2;
        }

        bool ValidAddGlcNAc()
        {
            return table_[0] == 2 && table_[4] > 0;
        }

        NGlycanHybrid CreateByAddGlcNAcCore()
        {
            var g = new NGlycanHybrid();
            g.SetTable(table_);
            g.table_[0] = g.table_[0] + 1;
            g.SetComposition(composite);
            g.AddMonosaccharide(Monosaccharide.GlcNAc);

            return g;
        }

        bool ValidAddGlcNAcBisect()
        {
            //bisect 0, not extanding on GlcNAc
            return table_[1] == 1 && table_[3] == 0 && table_[4] == 0;
        }

        NGlycanHybrid CreateByAddGlcNAcBisect()
        {
            var g = new NGlycanHybrid();
            g.SetTable(table_);
            g.table_[3] = 1;
            g.SetComposition(composite);
            g.AddMonosaccharide(Monosaccharide.GlcNAc);
            return g;
        }

        bool ValidAddGlcNAcBranch()
        {
            if (table_[5] == 0)
                return false;

            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
                {
                    if (table_[i + 8] == table_[i + 10] && table_[i + 12] == 0 && table_[i + 14] == 0 && table_[i + 16] == 0)
                    //equal GlcNAc Gal, no Fucose attached at terminal, no terminal NeuAc, NeuGc
                    {
                        return true;
                    }
                }
            }
            return false;

        }

        List<NGlycanHybrid> CreateByAddGlcNAcBranch()
        {
            List<NGlycanHybrid> glycans = new List<NGlycanHybrid>();
            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
                {
                    if (table_[i + 8] == table_[i + 10] && table_[i + 12] == 0 && table_[i + 14] == 0 && table_[i + 16] == 0)
                    {
                        var g = new NGlycanHybrid();
                        g.SetTable(table_);
                        g.table_[i + 8] = g.table_[i + 8] + 1;
                        g.SetComposition(composite);
                        g.AddMonosaccharide(Monosaccharide.GlcNAc);
                        glycans.Add(g);
                    }
                }
            }
            return glycans;
        }

        bool ValidAddManCore()
        {
            return table_[0] == 2 && (table_[4] == 0 || table_[5] == 0);
        }

        List<NGlycanHybrid> CreateByAddManCore()
        {
            List<NGlycanHybrid> glycans = new List<NGlycanHybrid>();
            if (table_[1] == 0)
            {
                var g = new NGlycanHybrid();
                g.SetTable(table_);
                g.SetComposition(composite);
                g.table_[1] = 1;
                g.AddMonosaccharide(Monosaccharide.Man);
                glycans.Add(g);
            }
            else
            {
                for(int i = 0; i < 2; i++)
                {
                    if (table_[4 + i] == 0)
                    {
                        var g = new NGlycanHybrid();
                        g.SetTable(table_);
                        g.SetComposition(composite);
                        g.table_[4 + i] = 1;
                        g.AddMonosaccharide(Monosaccharide.Man);
                        glycans.Add(g);
                    }
                }
            }
            
            return glycans;
        }

        bool ValidAddManBranch()
        {
            return table_[4] > 0;
        }

        List<NGlycanHybrid> CreateByAddManBranch()
        {
            List<NGlycanHybrid> glycans = new List<NGlycanHybrid>();
            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 6] < table_[i + 5]) // make it order
                {
                    var g = new NGlycanHybrid();
                    g.SetTable(table_);
                    g.table_[i + 6] = g.table_[i + 6] + 1;
                    g.SetComposition(composite);
                    g.AddMonosaccharide(Monosaccharide.Man);
                    glycans.Add(g);
                }
            }
            return glycans;
        }

        bool ValidAddGal()
        {
            if (table_[5] == 0)
                return false;

            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
                {
                    if (table_[i + 6] == table_[i + 8] + 1)
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        List<NGlycanHybrid> CreateByAddGal()
        {
            List<NGlycanHybrid> glycans = new List<NGlycanHybrid>();
            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
                {
                    if (table_[i + 6] == table_[i + 8] + 1)
                    {
                        var g = new NGlycanHybrid();
                        g.SetTable(table_);
                        g.table_[i + 8] = g.table_[i + 8] + 1;
                        g.SetComposition(composite);
                        g.AddMonosaccharide(Monosaccharide.Gal);
                        glycans.Add(g);
                    }
                }
            }
            return glycans;
        }

        bool ValidAddFucCore()
        {
            return table_[0] >= 1 && table_[2] == 0;  //core
        }

        NGlycanHybrid CreateByAddFucCore()
        {
            var g = new NGlycanHybrid();
            g.SetTable(table_);
            g.table_[2] = 1;
            g.SetComposition(composite);
            g.AddMonosaccharide(Monosaccharide.Fuc);
            return g;
        }

        bool ValidAddFucTerminal()
        {
            if (table_[5] == 0)
                return false;

            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 12] < table_[i + 11]) // make it order
                {
                    if (table_[i + 12] == 0 && table_[i + 8] > 0)
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        List<NGlycanHybrid> CreateByAddFucTerminal()
        {
            List<NGlycanHybrid> glycans = new List<NGlycanHybrid>();
            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 12] < table_[i + 11]) // make it order
                {
                    if (table_[i + 12] == 0 && table_[i + 8] > 0)
                    {
                        var g = new NGlycanHybrid();
                        g.SetTable(table_);
                        g.table_[i + 12] = 1;
                        g.SetComposition(composite);
                        g.AddMonosaccharide(Monosaccharide.Fuc);
                        glycans.Add(g);
                    }
                }
            }
            return glycans;
        }

        bool ValidAddNeuAc()
        {
            if (table_[5] == 0)
                return false;

            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 14] < table_[i + 13]) // make it order
                {
                    if (table_[i + 8] > 0 && table_[i + 8] == table_[i + 10] && table_[i + 14] == 0 && table_[i + 16] == 0)
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        List<NGlycanHybrid> CreateByAddNeuAc()
        {
            List<NGlycanHybrid> glycans = new List<NGlycanHybrid>();
            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 18] < table_[i + 13]) // make it order
                {
                    if (table_[i + 8] > 0 && table_[i + 8] == table_[i + 10] && table_[i + 14] == 0 && table_[i + 16] == 0)
                    {
                        var g = new NGlycanHybrid();
                        g.SetTable(table_);
                        g.table_[i + 14] = 1;
                        g.SetComposition(composite);
                        g.AddMonosaccharide(Monosaccharide.NeuAc);
                        glycans.Add(g);
                    }
                }
            }
            return glycans;
        }

        bool ValidAddNeuGc()
        {
            if (table_[5] == 0)
                return false;

            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 16] < table_[i + 15]) // make it order
                {
                    if (table_[i + 8] > 0 && table_[i + 8] == table_[i + 10] && table_[i + 14] == 0 && table_[i + 16] == 0)
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        List<NGlycanHybrid> CreateByAddNeuGc()
        {
            List<NGlycanHybrid> glycans = new List<NGlycanHybrid>();
            for (int i = 0; i < 2; i++)
            {
                if (i == 0 || table_[i + 16] < table_[i + 15]) // make it order
                {
                    if (table_[i + 8] > 0 && table_[i + 8] == table_[i + 10] && table_[i + 14] == 0 && table_[i + 16] == 0)
                    {
                        var g = new NGlycanHybrid();
                        g.SetTable(table_);
                        g.table_[i + 16] = 1;
                        g.SetComposition(composite);
                        g.AddMonosaccharide(Monosaccharide.NeuGc);
                        glycans.Add(g);
                    }
                }
            }
            return glycans;
        }

        public override GlycanType Type()
        {
            return GlycanType.NGlycanHybrid;
        }
    }
}
