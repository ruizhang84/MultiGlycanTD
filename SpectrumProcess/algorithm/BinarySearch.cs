using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpectrumProcess.algorithm
{
    public class BinarySearch<T> : ISearch<T>
    {
        protected double tolerance_;
        protected ToleranceBy type_;
        protected List<Point<T>> data_ = new List<Point<T>>();

        public BinarySearch(ToleranceBy by, double tol)
        {
            tolerance_ = tol;
            type_ = by;
        }

        public void Add(Point<T> point)
        {
            data_.Add(point);
        }

        public void Init(List<Point<T>> inputs)
        {
            data_ = inputs;
            data_.Sort();
        }

        public void Init()
        {
            data_.Sort();
        }

        public bool Match(double expect, double baseValue)
        {
            return BinarySearchPoints(expect, baseValue) >= 0;
        }

        public bool Match(double expect)
        {
            return Match(expect, expect);
        }

        public List<Point<T>> Search(double expect, double baseValue)
        {
            return ExtendAllMatch(expect, BinarySearchPoints(expect, baseValue), baseValue);
        }

        bool IsMatch(double expect, double observe, double baseValue)
        {
            if (type_ == ToleranceBy.PPM)
            {
                return Math.Abs(expect - observe) / baseValue * 1000000.0 < tolerance_;
            }
            return Math.Abs(expect - observe) < tolerance_;
        }

        public int BinarySearchPoints(double expect, double baseValue)
        {
            int start = 0;
            int end = data_.Count - 1;

            while (start <= end)
            {
                int mid = end + (start - end) / 2;
                if (IsMatch(data_[mid].Value(), expect, baseValue))
                {
                    return mid;
                }
                else if (data_[mid].Value() > expect)
                {
                    end = mid - 1;
                }
                else
                {
                    start = mid + 1;
                }
            }

            return -1;
        }

        public List<Point<T>> ExtendAllMatch(double expect, int matchIndx, double baseValue)
        {
            List<Point<T>> searched = new List<Point<T>>();
            if (matchIndx < 0) return searched;

            for (int left = matchIndx; left >= 0 && 
                IsMatch(expect, data_[left].Value(), baseValue); left--)
            {
                searched.Add(data_[left]);
            }

            for (int right = matchIndx + 1; right < data_.Count &&
                IsMatch(expect, data_[right].Value(), baseValue); right++)
            {
                searched.Add(data_[right]);
            }

            return searched;
        }

        public List<Point<T>> Search(double expect)
        {
            return Search(expect, expect);
        }

        public List<T> SearchContent(double expect, double baseValue)
        {
            return Search(expect, baseValue).Select(r => r.Content()).ToList();
        }

        public List<T> SearchContent(double expect)
        {
            return SearchContent(expect, expect);
        }

        public void SetTolerance(double tol)
        {
            tolerance_ = tol;
        }

        public void SetToleranceBy(ToleranceBy by)
        {
            type_ = by;
        }

        public double Tolerance()
        {
            return tolerance_;
        }

        public ToleranceBy ToleranceType()
        {
            return type_;
        }
    }
}
