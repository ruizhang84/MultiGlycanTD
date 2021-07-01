using MultiGlycanTDLibrary.model.glycan;
using MultiGlycanClassLibrary.util.mass;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using SpectrumProcess.brain;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public class GlycanBuilderFiltered : GlycanBuilder, IGlycanBuilder
    {
        List<SortedDictionary<Monosaccharide, int>> filterList;
        HashSet<string> filterSet;

        public GlycanBuilderFiltered(List<SortedDictionary<Monosaccharide, int>> filterList,
            int hexNAc = 12, int hex = 12, int fuc = 5, int neuAc = 4, int neuGc = 0,
            bool complex = true, bool hybrid = false, bool highMannose = false, int order = 10,
            bool permethylated = true, bool reduced = true, 
            Derivatization derivatization = Derivatization.Underivatized, int thread = 4) : 
            base(hexNAc, hex, fuc, neuAc, neuGc, complex, 
                hybrid, highMannose, order, 
                permethylated, reduced, derivatization, thread)
        {
            this.filterList = filterList;
            filterSet = new HashSet<string>();
            foreach (SortedDictionary<Monosaccharide, int> pairs in filterList)
            {
                filterSet.Add(CompositionToString(ConvertComposition(pairs)));
            }
        }

        public SortedDictionary<Monosaccharide, int> ConvertComposition(
            SortedDictionary<Monosaccharide, int> composition)
        {
            SortedDictionary<Monosaccharide, int> simple =
                new SortedDictionary<Monosaccharide, int>();
            foreach (Monosaccharide sugar in composition.Keys)
            {
                switch (sugar)
                {
                    case Monosaccharide.GlcNAc:
                        if (!simple.ContainsKey(Monosaccharide.HexNAc))
                            simple[Monosaccharide.HexNAc] = 0;
                        simple[Monosaccharide.HexNAc] += composition[Monosaccharide.GlcNAc];
                        break;
                    case Monosaccharide.Man:
                        if (!simple.ContainsKey(Monosaccharide.Hex))
                            simple[Monosaccharide.Hex] = 0;
                        simple[Monosaccharide.Hex] += composition[Monosaccharide.Man];
                        break;
                    case Monosaccharide.Gal:
                        if (!simple.ContainsKey(Monosaccharide.Hex))
                            simple[Monosaccharide.Hex] = 0;
                        simple[Monosaccharide.Hex] += composition[Monosaccharide.Gal];
                        break;
                    case Monosaccharide.HexNAc:
                        if (!simple.ContainsKey(Monosaccharide.HexNAc))
                            simple[Monosaccharide.HexNAc] = 0;
                        simple[Monosaccharide.HexNAc] += composition[Monosaccharide.HexNAc];
                        break;
                    case Monosaccharide.Hex:
                        if (!simple.ContainsKey(Monosaccharide.Hex))
                            simple[Monosaccharide.Hex] = 0;
                        simple[Monosaccharide.Hex] += composition[Monosaccharide.Hex];
                        break;
                    case Monosaccharide.Fuc:
                        simple[Monosaccharide.Fuc] = composition[Monosaccharide.Fuc];
                        break;
                    case Monosaccharide.NeuAc:
                        simple[Monosaccharide.NeuAc] = composition[Monosaccharide.NeuAc];
                        break;
                    case Monosaccharide.NeuGc:
                        simple[Monosaccharide.NeuGc] = composition[Monosaccharide.NeuGc];
                        break;
                    default:
                        break;
                }
            }
            return simple;
        }

        public string CompositionToString(
            SortedDictionary<Monosaccharide, int> composition)
        {
            string name = "";
            foreach (Monosaccharide sugar in composition.Keys)
            {
                switch (sugar)
                {
                    case Monosaccharide.HexNAc:
                        name += "HexNAc-" + composition[sugar] + " ";
                        break;
                    case Monosaccharide.Hex:
                        name += "Hex-" + composition[sugar] + " ";
                        break;
                    case Monosaccharide.Fuc:
                        name += "Fuc-" + composition[sugar] + " ";
                        break;
                    case Monosaccharide.NeuAc:
                        name += "NeuAc-" + composition[sugar] + " ";
                        break;
                    case Monosaccharide.NeuGc:
                        name += "NeuGc-" + composition[sugar] + " ";
                        break;
                    default:
                        break;
                }
            }
            return name;
        }

        public override Dictionary<string, IGlycan> GlycanMaps()
        {
            Dictionary<string, IGlycan> results = new Dictionary<string, IGlycan>();
            foreach (string id in glycans_map_.Keys)
            {
                IGlycan g = glycans_map_[id];
                if (g.IsValid() &&
                    filterSet.Contains(
                        CompositionToString(ConvertComposition(g.Composition()))))
                    results[id] = g;
            }
            return results;
        }

        protected override void BuildMaps()
        {
            foreach (var pair in glycans_map_)
            {
                IGlycan g = pair.Value;
                if (!g.IsValid() ||
                    !filterSet.Contains(
                        CompositionToString(ConvertComposition(g.Composition()))))
                    continue;
                string name = g.Name();
                if (!compound_map_.ContainsKey(name))
                {
                    compound_map_[name] = BuildCompound(g);
                    glycan_compound_map_[name] = new List<IGlycan>();
                }
                glycan_compound_map_[name].Add(g);
            };

            Parallel.ForEach(compound_map_, pair =>
            {
                string name = pair.Key;
                Compound compound = pair.Value;
                distr_map_[name] = Brain.Run.Distribute(compound, Order);
                mass_map_[name] = Brain.Run.CenterMass(compound, Order);
            });
        }

    }
}
