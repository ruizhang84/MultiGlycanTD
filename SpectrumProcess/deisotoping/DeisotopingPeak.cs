using SpectrumData;
using SpectrumData.Spectrum;

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

        public DeisotopingPeak(IPeak peak) :
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
