using MultiGlycanTDLibrary.model.glycan;
using MultiGlycanClassLibrary.util.mass;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MultiGlycanTDLibrary.util.brain;
using System.Collections.Concurrent;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public class GlycanBuilder
    {
        protected int hexNAc_;
        protected int hex_;
        protected int fuc_;
        protected int neuAc_;
        protected int neuGc_;
        public bool ComplexInclude { get; set; }
        public bool HybridInclude { get; set; }
        public bool HighMannoseInclude { get; set; }
        public bool Permethylated { get; set; } = true;
        public int Order { get; set; } = 10;

        public bool Reduced { get; set; } = true;

        protected Dictionary<string, IGlycan> glycans_map_; // glycan id -> glycan
        protected List<Monosaccharide> candidates_;

        protected Dictionary<string, List<IGlycan>> glycan_compound_map_;
        protected Dictionary<string, Compound> compound_map_;
        protected ConcurrentDictionary<string, List<double>> distr_map_;
        protected ConcurrentDictionary<string, List<double>> mass_map_;

        public GlycanBuilder(int hexNAc = 12, int hex = 12, int fuc = 5, int neuAc = 4, int neuGc = 0,
            bool complex = true, bool hybrid = false, bool highMannose = false, int order = 10,
            bool permethylated = true, bool reduced = true)
        {
            hexNAc_ = hexNAc;
            hex_ = hex;
            fuc_ = fuc;
            neuAc_ = neuAc;
            neuGc_ = neuGc;
            ComplexInclude = complex;
            HybridInclude = hybrid;
            HighMannoseInclude = highMannose;
            Permethylated = permethylated;
            Order = order;
            Reduced = reduced;
            compound_map_ = new Dictionary<string, Compound>();
            glycan_compound_map_ = new Dictionary<string, List<IGlycan>>();
            distr_map_ = new ConcurrentDictionary<string, List<double>>();
            mass_map_ = new ConcurrentDictionary<string, List<double>>();

            glycans_map_ = new Dictionary<string, IGlycan>();
            candidates_ = new List<Monosaccharide>()
            {
                Monosaccharide.GlcNAc,
                Monosaccharide.Man,
                Monosaccharide.Gal,
                Monosaccharide.Fuc,
                Monosaccharide.NeuAc,
                Monosaccharide.NeuGc
            };
        }

        public Dictionary<string, IGlycan> GlycanMaps()
        { return glycans_map_; }

        public Dictionary<string, List<IGlycan>> GlycanCompositionMaps()
        { return glycan_compound_map_;  }

        public Dictionary<string, List<double>> GlycanDistribMaps()
        { return distr_map_.ToDictionary(entry => entry.Key, entry => entry.Value); }

        public Dictionary<string, List<double>> GlycanMasMaps()
        { return mass_map_.ToDictionary(entry => entry.Key, entry => entry.Value);  }

        public void Build()
        {
            Queue<IGlycan> queue = new Queue<IGlycan>();
            IGlycan root;

            if (ComplexInclude)
            {
                root = new NGlycanComplex();
                queue.Enqueue(root);
            }

            if (HybridInclude)
            {
                root = new NGlycanHybrid();
                queue.Enqueue(root);
            }

            if (HighMannoseInclude)
            {
                root = new NHighMannose();
                queue.Enqueue(root);
            }

            while (queue.Count > 0)
            {
                IGlycan node = queue.Peek();
                queue.Dequeue();

                // next
                foreach (var it in candidates_)
                {
                    List<IGlycan> res = node.Grow(it);

                    foreach (var g in res)
                    {
                        if (SatisfyCriteria(g))
                        {
                            string id = g.ID();
                            if (!glycans_map_.ContainsKey(id))
                            {
                                g.Add(node);
                                glycans_map_[id] = g;
                                string name = g.Name();
                                if (!compound_map_.ContainsKey(name))
                                {
                                    compound_map_[name] = BuildCompound(g);
                                    glycan_compound_map_[name] = new List<IGlycan>();
                                }
                                glycan_compound_map_[name].Add(g);
                                queue.Enqueue(glycans_map_[id]);
                            }
                            else
                            {
                                glycans_map_[id].Add(node);
                            }
                        }
                    }
                }

            }

            Parallel.ForEach(compound_map_, pair =>
            {
                string name = pair.Key;
                Compound compound = pair.Value;
                distr_map_[name] = Brain.Run.Distribute(compound, Order);
                mass_map_[name] = Brain.Run.CenterMass(compound, Order);
            });


        }

        bool SatisfyCriteria(IGlycan glycan)
        {
            int hexNAc = 0, hex = 0, fuc = 0, neuAc = 0, neuGc = 0;
            SortedDictionary<Monosaccharide, int> composite = glycan.Composition();
            foreach (var key in composite.Keys)
            {
                switch (key)
                {
                    case Monosaccharide.GlcNAc:
                        hexNAc += composite[key];
                        break;
                    case Monosaccharide.Gal:
                        hex += composite[key];
                        break;
                    case Monosaccharide.Man:
                        hex += composite[key];
                        break;
                    case Monosaccharide.Fuc:
                        fuc += composite[key];
                        break;
                    case Monosaccharide.NeuAc:
                        neuAc += composite[key];
                        break;
                    case Monosaccharide.NeuGc:
                        neuGc += composite[key];
                        break;
                    default:
                        break;
                }
            }
            return (hexNAc <= hexNAc_ && hex <= hex_ && fuc <= fuc_
                    && neuAc <= neuAc_ && neuGc <= neuGc_);
        }

        public Compound BuildCompound(IGlycan glycan)
        {
            Dictionary<ElementType, int> formulaComposition = new Dictionary<ElementType, int>();
            var compose = glycan.Composition();

            // composition to elements
            foreach (var sugar in compose.Keys)
            {
                Dictionary<ElementType, int> tempCompose = NMonosaccharideCreator.Get.SubCompositions(
                    sugar, Permethylated);
                foreach (ElementType elm in tempCompose.Keys)
                {
                    if (!formulaComposition.ContainsKey(elm))
                    {
                        formulaComposition[elm] = 0;
                    }
                    formulaComposition[elm] += tempCompose[elm] * compose[sugar];
                }
            }

            // reducing end
            if (Permethylated)
            {
                if(Reduced)
                {
                    formulaComposition[ElementType.C] += 3;
                    formulaComposition[ElementType.H] += 10;
                    formulaComposition[ElementType.O] += 1;
                }
                else
                {
                    formulaComposition[ElementType.C] += 2;
                    formulaComposition[ElementType.H] += 6;
                    formulaComposition[ElementType.O] += 1;
                }
            }

            Compound formula = new Compound(formulaComposition);
            return formula;
        }
    }
}
