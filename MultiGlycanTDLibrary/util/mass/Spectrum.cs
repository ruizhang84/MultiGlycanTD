using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.util.mass
{
    public class Spectrum
    {
        protected static readonly Lazy<Spectrum>
            lazy = new Lazy<Spectrum>(() => new Spectrum());

        public static Spectrum To { get { return lazy.Value; } }

        protected List<double> ions;
        public static readonly double Proton = 1.0078;
        public static readonly double Ammonium = 14.00307 + 1.0078 * 4;
        public static readonly double Sodium = 22.98977;

        protected Spectrum()
        {
            ions = new List<double> { Proton, Ammonium, Sodium };
        }

        public void SetChargeIons(List<double> ionMass)
        {
            ions = ionMass;
        }

        public double ComputePPM(double expected, double observed)
        {
            return Math.Abs(expected - observed) / expected * 1000000.0;
        }

        public double Compute(double mz, double ion, int charge)
        {
            return (mz - ion) * charge;
        }

        public double ComputeMZ(double mass, double ion, int charge)
        {
            return (mass + ion * charge) / charge;
        }

    }
}
