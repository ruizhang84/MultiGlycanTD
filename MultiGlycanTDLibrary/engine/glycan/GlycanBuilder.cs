﻿using MultiGlycanTDLibrary.model.glycan;
using MultiGlycanClassLibrary.util.mass;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MultiGlycanTDLibrary.util.brain;

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

        protected Dictionary<string, IGlycan> glycans_map_; // glycan id -> glycan
        protected List<Monosaccharide> candidates_;
        protected bool permethylated = true;
        protected int order = 10;
        protected Dictionary<string, List<double>> mem_distr_;
        protected Dictionary<string, List<double>> mem_mass_;

        public GlycanBuilder(int hexNAc = 12, int hex = 12, int fuc = 5, int neuAc = 4, int neuGc = 0,
            bool complex = true, bool hybrid = false, bool highMannose = false, bool permethylated = true)
        {
            hexNAc_ = hexNAc;
            hex_ = hex;
            fuc_ = fuc;
            neuAc_ = neuAc;
            neuGc_ = neuGc;
            ComplexInclude = complex;
            HybridInclude = hybrid;
            HighMannoseInclude = highMannose;
            this.permethylated = permethylated;
            mem_distr_ = new Dictionary<string, List<double>>();
            mem_mass_ = new Dictionary<string, List<double>>();

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

        public List<Monosaccharide> Candidates() { return candidates_; }

        public int HexNAc() { return hexNAc_; }
        public int Hex() { return hex_; }
        public int Fuc() { return fuc_; }
        public int NeuAc() { return neuAc_; }
        public int NeuGc() { return neuGc_; }
        public void set_candidates(List<Monosaccharide> sugars)
        { candidates_ = sugars; }
        public void set_HexNAc(int num) { hexNAc_ = num; }
        public void set_Hex(int num) { hex_ = num; }
        public void set_Fuc(int num) { fuc_ = num; }
        public void set_NeuAc(int num) { neuAc_ = num; }
        public void set_NeuGc(int num) { neuGc_ = num; }

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
                                glycans_map_[id] = BuildDistribution(g);
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

        protected Compound BuildCompound(IGlycan glycan)
        {
            Dictionary<Element, int> formulaComposition = new Dictionary<Element, int>();
            var compose = glycan.Composition();
            foreach (var sugar in compose.Keys)
            {
                Dictionary<Element, int> tempCompose = NMonosaccharideCreator.Get.SubCompositions(
                    sugar, permethylated);
                foreach (Element elm in tempCompose.Keys)
                {
                    if (!formulaComposition.ContainsKey(elm))
                    {
                        formulaComposition[elm] = 0;
                    }
                    formulaComposition[elm] += tempCompose[elm] * compose[sugar];
                }
            }
            Compound formula = new Compound(formulaComposition);
            return formula;
        }

        protected IGlycan BuildDistribution(IGlycan glycan)
        {
            Compound formula = BuildCompound(glycan);
            glycan.SetFormula(formula);

            if (!mem_mass_.ContainsKey(formula.Name))
            {
                mem_distr_[formula.Name] = Brain.Run.Distribute(glycan.Formula(), order);
                mem_mass_[formula.Name] = Brain.Run.CenterMass(glycan.Formula(), order);
            }
            List<double> distrib = mem_distr_[formula.Name];
            List<double> massDistr = mem_mass_[formula.Name];

            glycan.SetDistrib(distrib);
            glycan.SetMass(massDistr);

            return glycan;
        }

    }
}
