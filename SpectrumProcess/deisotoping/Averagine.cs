using SpectrumProcess.brain;
using System;
using System.Collections.Generic;

namespace SpectrumProcess.deisotoping
{
    public enum AveragineType
    {
        Peptide,
        GlycoPeptide,
        Glycan,
        PermethylatedGlycan
    }

    public class Averagine
    {
        public AveragineType Type { get; set; }
        public Averagine(AveragineType type = AveragineType.PermethylatedGlycan)
        {
            Type = type;
        }

        Element CreateElement(ElementType type)
        {
            switch (type)
            {
                case ElementType.C:
                    return new C();
                case ElementType.H:
                    return new H();
                case ElementType.O:
                    return new O();
                case ElementType.N:
                    return new N();
                case ElementType.S:
                    return new S();
            }
            return new H();
        }

        public Compound Fit(double mass)
        {
            Dictionary<ElementType, int> compos
                = new Dictionary<ElementType, int>();

            // compute scale
            double total = 0;
            foreach (var item in Composition())
            {
                Element e = CreateElement(item.Key);
                double num = item.Value;
                total += e.Mass[e.Atomic] * num;
            }

            double scale = mass / total;

            // compute additional H
            total = 0;
            foreach (var item in Composition())
            {
                Element e = CreateElement(item.Key);
                int num = (int)Math.Floor(item.Value * scale);
                compos[item.Key] = num;
                total += e.Mass[e.Atomic] * num;
            }
            int addit = (int)Math.Floor(mass - total);
            foreach (var item in compos)
            {
                if (CreateElement(item.Key).Atomic == 1)
                {
                    compos[item.Key] = item.Value + addit;
                    break;
                }
            }

            return new Compound(compos);
        }

        public Dictionary<ElementType, double> Composition()
        {
            return Composition(Type);
        }

        public Dictionary<ElementType, double> Composition(AveragineType type)
        {
            switch (type)
            {
                case AveragineType.Peptide:
                    return Peptide;
                case AveragineType.GlycoPeptide:
                    return Glycopeptide;
                case AveragineType.Glycan:
                    return Glycan;
                case AveragineType.PermethylatedGlycan:
                    return PermethylatedGlycan;
                default:
                    break;
            }
            return Peptide;
        }

        public Dictionary<ElementType, double> Peptide
            = new Dictionary<ElementType, double>()
            {
                {ElementType.C, 4.9384 },
                {ElementType.H, 7.7583 },
                {ElementType.N, 1.3577 },
                {ElementType.O, 1.4773 },
                {ElementType.S, 0.0417 }
            };
        public Dictionary<ElementType, double> Glycopeptide
            = new Dictionary<ElementType, double>()
            {
                {ElementType.C, 10.93 },
                {ElementType.H, 15.75 },
                {ElementType.N, 1.6577 },
                {ElementType.O, 6.4773 },
                {ElementType.S, 0.02054 }
            };
        public Dictionary<ElementType, double> Glycan
            = new Dictionary<ElementType, double>()
            {
                {ElementType.C, 7.0 },
                {ElementType.H, 11.8333 },
                {ElementType.N, 0.5 },
                {ElementType.O, 5.16666 }
            };
        public Dictionary<ElementType, double> PermethylatedGlycan
            = new Dictionary<ElementType, double>()
            {
                {ElementType.C, 12.0 },
                {ElementType.H, 21.8333 },
                {ElementType.N, 0.5 },
                {ElementType.O, 5.16666 }
            };

    }
}
