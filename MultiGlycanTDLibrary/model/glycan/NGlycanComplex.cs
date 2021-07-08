using System.Collections.Generic;

namespace MultiGlycanTDLibrary.model.glycan
{
    public class NGlycanComplex : BaseNGlycan, IGlycan
    {
        //GlcNAc(2) - Man(1) - Fuc(1) - GlcNAc(bisect,1) -0,1,2,3
        //[Man(antenna1) - Man(antenna2)] - 4, 5,
        //[GlcNAc(branch1) - GlcNAc(branch2) - GlcNAc(branch3) - GlcNAc(branch4)] -6,7,8,9
        //[Gal(branch1) - Gal(branch2) - Gal(branch3) - Gal(branch4)] -10,11,12,13
        //[Fuc(branch1) - Fuc(branch2) - Fuc(branch3) - Fuc(branch4)] -14,15,16,17
        //[NeuAc(branch1) - NeuAc(branch2) - NeuAc(branch3) - NeuAc(branch4)] -18,19,20,21
        //[NeuGc(branch1) - NeuGc(branch2) - NeuGc(branch3) - NeuGc(branch4)] -22,23,24,25

        // 1) make it order based on observation that
        // a. fragments are not distinguish the bond linkage e.g. Y fragment
        // b. terminal position of fuc/neuAC/neuGc are ambiguous
        // 2) remove the sorting order if needed (recommend for small sugar).

        private int[] table_ = new int[26];

        public override IGlycan Clone()
        {
            IGlycan glycan = new NGlycanComplex();
            glycan.SetTable(table_);
            glycan.SetComposition(composite);
            return glycan;
        }

        public override bool IsValid()
        {
            // at least two chains
            if (table_[4] == 0 || table_[5] == 0)
                return false;
            if (table_[6] > 0 && table_[8] == 0)
                return false;
            // maksure sorted, so no duplicate
            if (table_[6] < table_[8])
                return false;
            if (table_[6] == table_[8] && table_[10] < table_[12])
                return false;
            if (table_[6] == table_[8] && table_[10] == table_[12]
                && (table_[14] < table_[16]
                || table_[18] < table_[20] || table_[22] < table_[24]))
                return false;
            for (int i = 0; i < 2; i++)
            {
                if (table_[i * 2 + 6] < table_[i * 2 + 7])
                    return false;
                if (table_[i * 2 + 6] == table_[i * 2 + 7]
                    && table_[i * 2 + 10] < table_[i * 2 + 11])
                    return false;
                if (table_[i * 2 + 6] == table_[i * 2 + 7]
                    && table_[i * 2 + 10] == table_[i * 2 + 11]
                    && (table_[i * 2 + 14] < table_[i * 2 + 15]
                    || table_[i * 2 + 18] < table_[i * 2 + 19]
                    || table_[i * 2 + 22] < table_[i * 2 + 23]))
                    return false;
            }

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
                        NGlycanComplex ptr = CreateByAddGlcNAcCore();
                        glycans.Add(ptr);
                    }
                    else if (ValidAddGlcNAc())
                    {
                        if (ValidAddGlcNAcBisect())
                        {
                            NGlycanComplex ptr = CreateByAddGlcNAcBisect();
                            glycans.Add(ptr);
                        }
                        if (ValidAddGlcNAcBranch())
                        {
                            List<NGlycanComplex> gs = CreateByAddGlcNAcBranch();
                            glycans.AddRange(gs);
                        }
                    }
                    break;

                case Monosaccharide.Man:
                    if (ValidAddMan())
                    {
                        List<NGlycanComplex> gs = CreateByAddMan();
                        glycans.AddRange(gs);
                    }
                    break;

                case Monosaccharide.Gal:
                    if (ValidAddGal())
                    {
                        List<NGlycanComplex> gs = CreateByAddGal();
                        glycans.AddRange(gs);
                    }
                    break;

                case Monosaccharide.Fuc:
                    if (ValidAddFucCore())
                    {
                        NGlycanComplex ptr = CreateByAddFucCore();
                        glycans.Add(ptr);
                    }
                    if (ValidAddFucTerminal())
                    {
                        List<NGlycanComplex> gs = CreateByAddFucTerminal();
                        glycans.AddRange(gs);
                    }
                    break;

                case Monosaccharide.NeuAc:
                    if (ValidAddNeuAc())
                    {
                        List<NGlycanComplex> gs = CreateByAddNeuAc();
                        glycans.AddRange(gs);
                    }
                    break;

                case Monosaccharide.NeuGc:
                    if (ValidAddNeuGc())
                    {
                        List<NGlycanComplex> gs = CreateByAddNeuGc();
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
            return table_[0] == 2 && (table_[4] == 1 || table_[5] == 1);
        }

        NGlycanComplex CreateByAddGlcNAcCore()
        {
            var g = new NGlycanComplex();
            g.SetTable(table_);
            g.table_[0] = g.table_[0] + 1;
            g.SetComposition(composite);
            g.AddMonosaccharide(Monosaccharide.GlcNAc);
            return g;
        }

        bool ValidAddGlcNAcBisect()
        {
            //bisect 0, not extanding on GlcNAc
            return table_[1] == 1 && table_[3] == 0;
        }

        NGlycanComplex CreateByAddGlcNAcBisect()
        {
            var g = new NGlycanComplex();
            g.SetTable(table_);
            g.table_[3] = 1;
            g.SetComposition(composite);
            g.AddMonosaccharide(Monosaccharide.GlcNAc);
            return g;
        }

        bool ValidAddGlcNAcBranch()
        {
            for (int i = 0; i < 2; i++)
            {
                // atenna attachable
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 6] == table_[i * 2 + j + 10] && table_[i * 2 + j + 14] == 0 && table_[i * 2 + j + 18] == 0 && table_[i * 2 + j + 22] == 0)
                        //equal GlcNAc Gal, no Fucose attached at terminal, no terminal NeuAc, NeuGc
                        {
                            return true;
                        }
                    }
                }
            }
            return false;

        }

        List<NGlycanComplex> CreateByAddGlcNAcBranch()
        {
            List<NGlycanComplex> glycans = new List<NGlycanComplex>();
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                // atenna attachable
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 6] == table_[i * 2 + j + 10] && table_[i * 2 + j + 14] == 0 && table_[i * 2 + j + 18] == 0 && table_[i * 2 + j + 22] == 0)
                        //equal GlcNAc Gal, no Fucose attached at terminal, no terminal NeuAc, NeuGc
                        {
                            var g = new NGlycanComplex();
                            g.SetTable(table_);
                            g.table_[i * 2 + j + 6] = g.table_[i * 2 + j + 6] + 1;
                            g.SetComposition(composite);
                            g.AddMonosaccharide(Monosaccharide.GlcNAc);
                            glycans.Add(g);
                        }
                    }
                }
            }
            return glycans;
        }

        bool ValidAddMan()
        {
            return table_[0] == 2 && (table_[4] == 0 || table_[5] == 0);
        }

        List<NGlycanComplex> CreateByAddMan()
        {
            List<NGlycanComplex> glycans = new List<NGlycanComplex>();
            if (table_[1] == 0)
            {
                var g = new NGlycanComplex();
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
                    if (table_[4 + i] == 0)
                    {
                        var g = new NGlycanComplex();
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

        bool ValidAddGal()
        {
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 6] == table_[i * 2 + j + 10] + 1)
                        {
                            return true;
                        }
                    }
                }
            }

            return false;
        }

        List<NGlycanComplex> CreateByAddGal()
        {
            List<NGlycanComplex> glycans = new List<NGlycanComplex>();
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 6] == table_[i * 2 + j + 10] + 1)
                        {
                            var g = new NGlycanComplex();
                            g.SetTable(table_);
                            g.table_[i * 2 + j + 10] = g.table_[i * 2 + j + 10] + 1;
                            g.SetComposition(composite);
                            g.AddMonosaccharide(Monosaccharide.Gal);
                            glycans.Add(g);
                        }
                    }
                }
            }
            return glycans;
        }

        bool ValidAddFucCore()
        {
            return table_[0] >= 1 && table_[2] == 0;  //core
        }

        NGlycanComplex CreateByAddFucCore()
        {
            var g = new NGlycanComplex();
            g.SetTable(table_);
            g.table_[2] = 1;
            g.SetComposition(composite);
            g.AddMonosaccharide(Monosaccharide.Fuc);
            return g;
        }

        bool ValidAddFucTerminal()
        {
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 14] == 0 && table_[i * 2 + j + 6] > 0)
                        {
                            return true;
                        }
                    }
                }
            }

            return false;
        }

        List<NGlycanComplex> CreateByAddFucTerminal()
        {
            List<NGlycanComplex> glycans = new List<NGlycanComplex>();
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 14] == 0 && table_[i * 2 + j + 6] > 0)
                        {
                            var g = new NGlycanComplex();
                            g.SetTable(table_);
                            g.table_[i * 2 + j + 14] = 1;
                            g.SetComposition(composite);
                            g.AddMonosaccharide(Monosaccharide.Fuc);
                            glycans.Add(g);
                        }
                    }
                }
            }

            return glycans;
        }

        bool ValidAddNeuAc()
        {
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 6] > 0 && table_[i * 2 + j + 6] == table_[i * 2 + j + 10] && table_[i * 2 + j + 18] == 0 && table_[i * 2 + j + 22] == 0)
                        {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        List<NGlycanComplex> CreateByAddNeuAc()
        {
            List<NGlycanComplex> glycans = new List<NGlycanComplex>();
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 6] > 0 && table_[i * 2 + j + 6] == table_[i * 2 + j + 10] && table_[i * 2 + j + 18] == 0 && table_[i * 2 + j + 22] == 0)
                        {
                            var g = new NGlycanComplex();
                            g.SetTable(table_);
                            g.table_[i * 2 + j + 18] = 1;
                            g.SetComposition(composite);
                            g.AddMonosaccharide(Monosaccharide.NeuAc);
                            glycans.Add(g);
                        }
                    }
                }
            }
            return glycans;
        }

        bool ValidAddNeuGc()
        {
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        if (table_[i * 2 + j + 6] > 0 && table_[i * 2 + j + 6] == table_[i * 2 + j + 10] && table_[i * 2 + j + 18] == 0 && table_[i * 2 + j + 22] == 0)
                        {
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        List<NGlycanComplex> CreateByAddNeuGc()
        {
            List<NGlycanComplex> glycans = new List<NGlycanComplex>();
            for (int i = 0; i < 2; i++)
            {
                if (table_[i + 4] > 0)
                {
                    for (int j = 0; j < 2; j++)
                    {

                        if (table_[i * 2 + j + 6] > 0 && table_[i * 2 + j + 6] == table_[i * 2 + j + 10] && table_[i * 2 + j + 18] == 0 && table_[i * 2 + j + 22] == 0)
                        {
                            var g = new NGlycanComplex();
                            g.SetTable(table_);
                            g.table_[i + 22] = 1;
                            g.SetComposition(composite);
                            g.AddMonosaccharide(Monosaccharide.NeuGc);
                            glycans.Add(g);
                        }
                    }
                }
            }
            return glycans;
        }

        public override GlycanType Type()
        {
            return GlycanType.NGlycanComplex;
        }
    }
}
