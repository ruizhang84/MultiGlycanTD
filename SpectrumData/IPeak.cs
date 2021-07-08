using System;

namespace SpectrumData
{
    public interface IPeak : IComparable<IPeak>
    {
        IPeak Clone();
        double GetMZ();
        void SetMZ(double mz);
        double GetIntensity();
        void SetIntensity(double intensity);
    }
}
