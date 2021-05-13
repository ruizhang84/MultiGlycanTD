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
        double Mass();
        void SetMass(double mass);
        IGlycan Clone();
        GlycanType Type();
        List<IGlycan> Children();
        void Add(IGlycan glycan);
        string Name();
        string ID();
        int[] Table();
        void SetTable(int[] table);
        SortedDictionary<Monosaccharide, int> Composition();
        void SetComposition(SortedDictionary<Monosaccharide, int> composition);
        List<IGlycan> Grow(Monosaccharide monosaccharide);
        bool IsValid();
        Compound Formula();
        void SetFormula(Compound formula);
        List<double> GetDistrib();
        void SetDistrib(List<double> distrib);
        // mass of glycan corresponding to most aboundant peak
        double HighestPeak();
        void SetHighestPeak(double peak);
    }
}
