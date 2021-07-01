using SpectrumData;
using SpectrumData.Spectrum;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpectrumProcess.deisotoping
{
    public class DeisotopingPeak : GeneralPeak, IPeak
    {
        protected int charge;
        protected bool hasCharge;
        public int Charge
        {
            get => charge;
            set
            {
                charge = value;
                hasCharge = true;
            }
        }

        public DeisotopingPeak(IPeak peak):
            base(peak.GetMZ(), peak.GetIntensity())
        {
            charge = 0;
            hasCharge = false;
        }

        public bool ChargeAssigned()
        {
            return hasCharge;
        }
    }
}
