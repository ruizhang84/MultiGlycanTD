using SpectrumData;
using System.Collections.Generic;

namespace SpectrumProcess
{
    public interface ICharger
    {
        int Charge(List<IPeak> peaks, double lower, double upper);
    }
}
