using MultiGlycanTDLibrary.model.glycan;
using SpectrumProcess.brain;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public enum Derivatization
    {
        Underivatized, k2AA, k2AB
    }

    public class GlycanBuilder : IGlycanBuilder
    {
        protected int hexNAc_;
        protected int hex_;
        protected int fuc_;
        protected int neuAc_;
        protected int neuGc_;
        public int Thread { get; set; }
        public bool ComplexInclude { get; set; }
        public bool HybridInclude { get; set; }
        public bool HighMannoseInclude { get; set; }
        public bool Permethylated { get; set; } = true;
        public bool Reduced { get; set; } = true;
        public int Order { get; set; } = 10;
        public Derivatization Derivates { get; set; }
            = Derivatization.Underivatized;

        protected Dictionary<string, IGlycan> glycans_map_; // glycan id -> glycan
        protected List<Monosaccharide> candidates_;

        protected Dictionary<string, List<IGlycan>> glycan_compound_map_;
        protected Dictionary<string, Compound> compound_map_;
        protected ConcurrentDictionary<string, List<double>> distr_map_;
        protected ConcurrentDictionary<string, List<double>> mass_map_;
        protected object obj = new object();

        public GlycanBuilder(int hexNAc = 12, int hex = 12, int fuc = 5, int neuAc = 4, int neuGc = 0,
            bool complex = true, bool hybrid = false, bool highMannose = false, int order = 10,
            bool permethylated = true, bool reduced = true,
            Derivatization derivatization = Derivatization.Underivatized, int thread = 4)
        {
            hexNAc_ = hexNAc;
            hex_ = hex;
            fuc_ = fuc;
            neuAc_ = neuAc;
            neuGc_ = neuGc;
            Thread = thread;
            ComplexInclude = complex;
            HybridInclude = hybrid;
            HighMannoseInclude = highMannose;
            Permethylated = permethylated;
            Reduced = reduced;
            Order = order;
            Derivates = derivatization;
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

        public virtual Dictionary<string, IGlycan> GlycanMaps()
        {
            Dictionary<string, IGlycan> results = new Dictionary<string, IGlycan>();
            foreach (string id in glycans_map_.Keys)
            {
                IGlycan g = glycans_map_[id];
                if (g.IsValid())
                    results[id] = g;
            }
            return results;
        }

        public Dictionary<string, List<IGlycan>> GlycanCompositionMaps()
        { return glycan_compound_map_; }

        public Dictionary<string, List<double>> GlycanDistribMaps()
        { return distr_map_.ToDictionary(entry => entry.Key, entry => entry.Value); }

        public Dictionary<string, List<double>> GlycanMassMaps()
        { return mass_map_.ToDictionary(entry => entry.Key, entry => entry.Value); }

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
                Queue<IGlycan> bags = new Queue<IGlycan>();
                Parallel.ForEach(queue, new ParallelOptions { MaxDegreeOfParallelism = Thread }, node =>
                {
                    // next
                    foreach (var it in candidates_)
                    {
                        List<IGlycan> res = node.Grow(it);
                        foreach (var g in res)
                        {
                            if (SatisfyCriteria(g))
                            {
                                string id = g.ID();
                                lock (obj)
                                {
                                    if (!glycans_map_.ContainsKey(id))
                                    {
                                        glycans_map_[id] = g;
                                        bags.Enqueue(glycans_map_[id]);
                                    }
                                }
                                glycans_map_[id].Add(node);
                            }
                        }
                    }
                });

                queue = bags;
            }

            BuildMaps();

        }

        protected virtual void BuildMaps()
        {
            foreach (var pair in glycans_map_)
            {
                IGlycan g = pair.Value;
                if (!g.IsValid())
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

        public bool SatisfyCriteria(IGlycan glycan)
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
                Dictionary<ElementType, int> tempCompose = NMonosaccharideCreator.Get.Compositions(
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
                if (Reduced)
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
            else
            // derivation
            {
                switch (Derivates)
                {
                    case Derivatization.k2AA:
                        formulaComposition[ElementType.C] += 7;
                        formulaComposition[ElementType.H] += 7;
                        formulaComposition[ElementType.O] += 1;
                        formulaComposition[ElementType.N] += 2;
                        break;
                    case Derivatization.k2AB:
                        formulaComposition[ElementType.C] += 7;
                        formulaComposition[ElementType.H] += 8;
                        formulaComposition[ElementType.O] += 2;
                        formulaComposition[ElementType.N] += 1;
                        break;
                    default:
                        formulaComposition[ElementType.H] += 2;
                        formulaComposition[ElementType.O] += 1;
                        break;
                }
            }

            Compound formula = new Compound(formulaComposition);
            return formula;
        }
    }
}
