using MultiGlycanTDLibrary.util.brain;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model.glycan
{
    public class NMonosaccharideCreator
    {
        protected static readonly Lazy<NMonosaccharideCreator>
            lazy = new Lazy<NMonosaccharideCreator>(() => new NMonosaccharideCreator());
        public static NMonosaccharideCreator Get { get { return lazy.Value; } }

        protected NMonosaccharideCreator() { }

        // composition in a glycan
        public Dictionary<ElementType, int> Compositions(Monosaccharide sugar, bool permethylated)
        {
            // GlcNAc C8H15NO6, Man/Gal C6H12O6, Fuc C6H12O5,  NeuAc C11H19NO9, NeuGc C11H19NO10
            if (permethylated)
            {
                switch (sugar)
                {
                    case Monosaccharide.GlcNAc:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 11}, {ElementType.H, 19}, {ElementType.N, 1}, {ElementType.O, 5} };
                    case Monosaccharide.Man:
                    case Monosaccharide.Gal:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 9}, {ElementType.H, 16}, {ElementType.O, 5} };
                    case Monosaccharide.Fuc:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 8}, {ElementType.H, 14}, {ElementType.O, 4} };
                    case Monosaccharide.NeuAc:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 16}, {ElementType.H, 27}, {ElementType.N, 1}, {ElementType.O, 8} };
                    case Monosaccharide.NeuGc:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 17}, {ElementType.H, 29}, {ElementType.N, 1}, {ElementType.O, 9} };
                }
            }
            else
            {
                switch (sugar)
                {
                    case Monosaccharide.GlcNAc:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 8}, {ElementType.H, 13}, {ElementType.N, 1}, {ElementType.O, 5} };
                    case Monosaccharide.Man:
                    case Monosaccharide.Gal:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 6}, {ElementType.H, 10}, {ElementType.O, 5} };
                    case Monosaccharide.Fuc:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 6}, {ElementType.H, 10}, {ElementType.O, 4} };
                    case Monosaccharide.NeuAc:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 11}, {ElementType.H, 17}, {ElementType.N, 1}, {ElementType.O, 8} };
                    case Monosaccharide.NeuGc:
                        return new Dictionary<ElementType, int>()
                    { {ElementType.C, 11}, {ElementType.H, 17}, {ElementType.N, 1}, {ElementType.O, 9} };
                }
            }



            return new Dictionary<ElementType, int>();
        }

    }
}
