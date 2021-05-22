using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace MultiGlycanTDLibrary.util.brain
{
    public class Brain
    {
        protected static readonly Lazy<Brain>
            lazy = new Lazy<Brain>(() => new Brain());
        public static Brain Run { get { return lazy.Value; } }
        protected Brain(){}

        public Element ElementConv(ElementType type)
        {
            switch(type)
            {
                case ElementType.C:
                    return new C();
                case ElementType.H:
                    return new H();
                case ElementType.N:
                    return new N();
                case ElementType.O:
                    return new O();
                case ElementType.S:
                    return new S();
            }
            return new H();
        }

        private double PsiFunc(Compound compound, int order)
        {
            Complex psi = 0;
            foreach(var item in compound.Composition)
            {
                int num = item.Value;
                Element element = ElementConv(item.Key);
                foreach (Complex root in element.Root())
                {
                    psi += num * Complex.Pow(root, -order);
                }
            }
            return psi.Real;
        }

        private double IteractionCoeff(List<double> coeff, List<double> polynom)
        {
            int j = coeff.Count;
            double qi = 0;
            for(int l = 0; l < j; l++)
            {
                qi += coeff[j - 1 - l] * polynom[l];
            }
            return -qi / j;
        }

        public List<double> Distribute(Compound compound, int order)
        {
            List<double> coeff = new List<double>();
            double q0 = 1.0;
            foreach(var item in compound.Composition)
            {
                int num = item.Value;
                Element element = ElementConv(item.Key);
                q0 *= Math.Pow(element.Abundance[element.Abundance.Keys.Min()], num);
            }
            coeff.Add(q0);


            double qi = q0;
            List<double> polynom = new List<double>();
            for(int i = 1; i < order; i++)
            {
                double psi = PsiFunc(compound, i);
                polynom.Add(psi);

                qi = IteractionCoeff(coeff, polynom);
                coeff.Add(qi);
            }
            return coeff;
        }

       public List<double> CenterMass(Compound compound, int order)
       {
            double[] variant = new double[order];
            Dictionary<ElementType, int> formular = compound.Composition;

            foreach (ElementType element in formular.Keys)
            {
                int count = formular[element];
                formular[element] -= 1;
                Element e = ElementConv(element);
                List<double> coeff = Distribute(new Compound(formular), order);

                switch (element)
                {
                    case ElementType.C:
                        for (int i = 0; i < coeff.Count; i++)
                        {
                            variant[i] += coeff[i] * count * e.Abundance[12] * e.Mass[12];
                            if (i < coeff.Count - 1)
                                variant[i + 1] += coeff[i] * count * e.Abundance[13] * e.Mass[13];
                        }
                        break;
                    case ElementType.H:
                        for (int i = 0; i < coeff.Count; i++)
                        {
                            variant[i] += coeff[i] * count * e.Abundance[1] * e.Mass[1];
                            if (i < coeff.Count - 1)
                                variant[i + 1] += coeff[i] * count * e.Abundance[2] * e.Mass[2];
                        }
                        break;
                    case ElementType.O:
                        for (int i = 0; i < coeff.Count; i++)
                        {
                            variant[i] += coeff[i] * count * e.Abundance[16] * e.Mass[16];
                            if (i < coeff.Count - 1)
                                variant[i + 1] += coeff[i] * count * e.Abundance[17] * e.Mass[17];
                            if (i < coeff.Count - 2)
                                variant[i + 2] += coeff[i] * count * e.Abundance[18] * e.Mass[18];
                        }
                        break;
                    case ElementType.N:
                        for (int i = 0; i < coeff.Count; i++)
                        {
                            variant[i] += coeff[i] * count * e.Abundance[14] * e.Mass[14];
                            if (i < coeff.Count - 1)
                                variant[i + 1] += coeff[i] * count * e.Abundance[15] * e.Mass[15];
                        }
                        break;
                    case ElementType.S:
                        for (int i = 0; i < coeff.Count; i++)
                        {
                            variant[i] += coeff[i] * count * e.Abundance[32] * e.Mass[32];
                            if (i < coeff.Count - 1)
                                variant[i + 1] += coeff[i] * count * e.Abundance[33] * e.Mass[33];
                            if (i < coeff.Count - 2)
                                variant[i + 2] += coeff[i] * count * e.Abundance[34] * e.Mass[34];
                            if (i < coeff.Count - 4)
                                variant[i + 4] += coeff[i] * count * e.Abundance[36] * e.Mass[36];
                        }
                        break;
                }
                formular[element] += 1;
            }
            List<double> coefficient = Distribute(new Compound(formular), order);
            for(int i = 0; i < order; i ++)
            {
                variant[i] /= coefficient[i];
            }

            return variant.ToList();
       }


    }
}
