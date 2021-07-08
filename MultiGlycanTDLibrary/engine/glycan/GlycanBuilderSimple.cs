using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.glycan
{
    public class GlycanBuilderSimple : IGlycanBuilder
    {
        protected int hexNAc_;
        protected int hex_;
        protected int fuc_;
        protected int neuAc_;
        protected int neuGc_;
        public int Thread { get; set; } = 4;
        public bool ComplexInclude { get; set; }
        public bool HybridInclude { get; set; }
        public bool HighMannoseInclude { get; set; }

        protected ConcurrentDictionary<string, IGlycan> glycans_map_; // glycan id -> glycan
        protected List<Monosaccharide> candidates_;

        protected object obj = new object();

        public GlycanBuilderSimple(int hexNAc = 12, int hex = 12, int fuc = 5, int neuAc = 4, int neuGc = 0,
            bool complex = true, bool hybrid = false, bool highMannose = false)
        {
            hexNAc_ = hexNAc;
            hex_ = hex;
            fuc_ = fuc;
            neuAc_ = neuAc;
            neuGc_ = neuGc;
            ComplexInclude = complex;
            HybridInclude = hybrid;
            HighMannoseInclude = highMannose;

            glycans_map_ = new ConcurrentDictionary<string, IGlycan>();
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

        }

        public Dictionary<string, List<IGlycan>> GlycanCompositionMaps()
        {
            throw new NotImplementedException();
        }

        public Dictionary<string, List<double>> GlycanDistribMaps()
        {
            throw new NotImplementedException();
        }

        public Dictionary<string, IGlycan> GlycanMaps()
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

        public Dictionary<string, List<double>> GlycanMassMaps()
        {
            throw new NotImplementedException();
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
    }
}
