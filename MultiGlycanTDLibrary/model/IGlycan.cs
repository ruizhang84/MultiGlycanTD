using MultiGlycanTDLibrary.util.brain;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.model.glycan
{
    public enum Monosaccharide
    { GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc }

    public enum GlycanType
    {
        NGlycanComplex, NGlycanHybrid, NHighMannose
    }

    public interface IGlycan
    {
        IGlycan Clone();
        bool IsValid();
        bool Sorted();
        GlycanType Type();
        HashSet<IGlycan> ChildrenHashSet();
        List<IGlycan> Children();
        List<IGlycan> Fragments();
        void Add(IGlycan glycan);
        string Name();
        string ID();
        int[] Table();
        void SetTable(int[] table);
        void SetSorted(bool sorted);
        SortedDictionary<Monosaccharide, int> Composition();
        void SetComposition(SortedDictionary<Monosaccharide, int> composition);
        List<IGlycan> Grow(Monosaccharide monosaccharide);

    }
}
