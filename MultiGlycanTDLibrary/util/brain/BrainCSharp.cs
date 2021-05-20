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

        private double PsiFunc(Compound compound, int order)
        {
            Complex psi = 0;
            foreach(var item in compound.Composition)
            {
                int num = item.Value;
                Element element = item.Key;
                foreach(Complex root in element.Root())
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
                Element element = item.Key;
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
            Dictionary<Element, int> formular = compound.Composition;
            Element C = new C();
            Element H = new H();
            Element N = new N();
            Element O = new O();
            Element S = new S();

            foreach(Element element in formular.Keys)
            {
                int count = formular[element];
                if (element.Name == "Carbon")
                {
                    formular[element] -= 1;
                    List<double> coeff = Distribute(new Compound(formular), order);
                    for(int i = 0; i < coeff.Count; i++)
                    {
                        variant[i] += coeff[i] * count * C.Abundance[12] * C.Mass[12];
                        if (i < coeff.Count - 1)
                            variant[i + 1] += coeff[i] * count * C.Abundance[13] * C.Mass[13];
                    }
                    formular[element] += 1;
                }
                else if (element.Name == "Hydrogen")
                {
                    formular[element] -= 1;
                    List<double> coeff = Distribute(new Compound(formular), order);
                    for (int i = 0; i < coeff.Count; i++)
                    {
                        variant[i] += coeff[i] * count * H.Abundance[1] * H.Mass[1];
                        if (i < coeff.Count - 1)
                            variant[i + 1] += coeff[i] * count * H.Abundance[2] * H.Mass[2];
                    }
                    formular[element] += 1;
                }
                else if (element.Name == "Nitrogen")
                {
                    formular[element] -= 1;
                    List<double> coeff = Distribute(new Compound(formular), order);
                    for (int i = 0; i < coeff.Count; i++)
                    {
                        variant[i] += coeff[i] * count * N.Abundance[14] * N.Mass[14];
                        if (i < coeff.Count - 1)
                            variant[i + 1] += coeff[i] * count * N.Abundance[15] * N.Mass[15];
                    }
                    formular[element] += 1;
                }
                else if (element.Name == "Oxygen")
                {
                    formular[element] -= 1;
                    List<double> coeff = Distribute(new Compound(formular), order);
                    for (int i = 0; i < coeff.Count; i++)
                    {
                        variant[i] += coeff[i] * count * O.Abundance[16] * O.Mass[16];
                        if (i < coeff.Count - 1)
                            variant[i + 1] += coeff[i] * count * O.Abundance[17] * O.Mass[17];
                        if (i < coeff.Count - 2)
                            variant[i + 2] += coeff[i] * count * O.Abundance[18] * O.Mass[18];
                    }
                    formular[element] += 1;
                }
                else if (element.Name == "Sulfur")
                {
                    formular[element] -= 1;
                    List<double> coeff = Distribute(new Compound(formular), order);
                    for (int i = 0; i < coeff.Count; i++)
                    {
                        variant[i] += coeff[i] * count * S.Abundance[32] * S.Mass[32];
                        if (i < coeff.Count - 1)
                            variant[i + 1] += coeff[i] * count * S.Abundance[33] * S.Mass[33];
                        if (i < coeff.Count - 2)
                            variant[i + 2] += coeff[i] * count * S.Abundance[34] * S.Mass[34];
                        if (i < coeff.Count - 4)
                            variant[i + 4] += coeff[i] * count * S.Abundance[36] * S.Mass[36];
                    }
                    formular[element] += 1;
                }
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
