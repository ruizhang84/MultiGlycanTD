using System.Collections.Generic;

namespace MultiGlycanTDLibrary.model.glycan
{
    public enum Monosaccharide
    { GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc, HexNAc, Hex }

    public enum GlycanType
    {
        NGlycanComplex, NGlycanHybrid, NHighMannose
    }

    public interface IGlycan
    {
        IGlycan Clone();
        bool IsValid(); // distinguish glycan and fragments, valid for glycan (i.e., pentacore)
        GlycanType Type();
        List<IGlycan> Children();
        List<IGlycan> Fragments();
        void Add(IGlycan glycan);
        string Name();
        string ID();
        int[] Table();
        void SetTable(int[] table);
        SortedDictionary<Monosaccharide, int> Composition();
        void SetComposition(SortedDictionary<Monosaccharide, int> composition);
        List<IGlycan> Grow(Monosaccharide monosaccharide);

    }
}
