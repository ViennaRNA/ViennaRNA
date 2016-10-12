#ifndef _ALIGNMENT_TABLES_H_
#define _ALIGNMENT_TABLES_H_

#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>

// superclass of tables, has the row start info

template<class R> 
class TAD_DP_Table {
	public: 
		friend std::ostream& operator<<(std::ostream &out, const TAD_DP_Table<R> &table) {
			table.print(out);
			return out;
		}

		TAD_DP_Table(unsigned long rows, unsigned long cols, R init) 
			: rows_(rows),
			cols_(cols),
			mtrxSize_(rows*cols) {
	    rowStart_ = new unsigned long[rows];
	    rowStart_[0] = 0;
	    for (unsigned long h = 1; h < rows; h++) {
	        rowStart_[h] = rowStart_[h - 1] + cols;
	    }
			//TODO if (topdown)
			computed_ = new bool[mtrxSize_];
		}

		~TAD_DP_Table(){
			delete[] rowStart_;
		}

		virtual void checkSpaceConsumption() = 0;
    virtual void print(std::ostream &s) const = 0;

		// TODO if nicht topdown dann was?
	  inline bool computed(const unsigned long i, const unsigned long j) const {
        assert(this->rowStart_[i] + j < this->mtrxSize_);
        return computed_[this->rowStart_[i] + j];
    };

    inline void setComputed(const unsigned long i, const unsigned long j) {
      assert(this->rowStart_[i] + j < this->mtrxSize_);
      computed_[this->rowStart_[i] + j] = true;
    };


	protected:
		unsigned long rows_;
		unsigned long cols_;
    unsigned long mtrxSize_;
    unsigned long *rowStart_;
		bool *computed_;

};


// table for linear alignment

template<class R> 
class TAD_DP_TableLinear : public TAD_DP_Table<R> {
	private:
    R *mtrx_;

	public:
		TAD_DP_TableLinear(unsigned long rows, unsigned long cols, R init) 
			: TAD_DP_Table<R>(rows,cols,init) {
			checkSpaceConsumption();
			mtrx_ = new R[this->mtrxSize_];
		}

		~TAD_DP_TableLinear() {
			delete[] mtrx_;
		}

		void checkSpaceConsumption() {
		    // check for an overflow
		    if (this->rows_ > ULONG_MAX / this->cols_) {
		        std::cerr << "Error: Overflow in calculation matrix multiplication. Calculation terminated." << std::endl;
		        exit(EXIT_FAILURE);
		    }
		    // maximum array size is 2GB
		    if (this->mtrxSize_ > 2000000000 || this->rows_ > 2000000000) {
		        std::cerr << "Error: Maximum array size of 2GB exceeded due to large input data. Calculation terminated." << std::endl;
		        exit(EXIT_FAILURE);
		    }
		}

    inline R getMtrxVal(const unsigned long i, const unsigned long j) const {
        assert(this->rowStart_[i] + j < this->mtrxSize_);
        return mtrx_[this->rowStart_[i] + j];
		}

		inline void setMtrxVal(const unsigned long i, const unsigned long j, R& val) {
      assert(this->rowStart_[i] + j < this->mtrxSize_);
      mtrx_[this->rowStart_[i] + j] = val;
		}

    void print(std::ostream &s) const {
			for (unsigned int i = 0; i < this->rows_; i++) {
				for (unsigned int j = 0; j < this->cols_; j++) {
					 s << mtrx_[this->rowStart_[i] + j] << " ";
				}
				s << std::endl;
			}
			s << std::endl;
		}
};


// tables for affine alignment

const int S = 0, V = 1, H = 2, V_ = 3, H_ = 4, V_H = 5, VH_ = 6;
const std::string table_name[] =  {"S","V","H","V'","H'","V'H","VH'"};

template<class R> 
class TAD_DP_TableAffine : public TAD_DP_Table<R> {
	private:
    R *mtrxS_;
    R *mtrxV_;
    R *mtrxH_;
    R *mtrxV__;
    R *mtrxH__;
    R *mtrxV_H_;
    R *mtrxVH__;
		int localOptimumTable_;

    inline R* getMtrx(int table) const {
        switch (table) {
          case S:
            return mtrxS_;
          case V:
            return mtrxV_;
          case H:
            return mtrxH_;
          case V_:
            return mtrxV__;
          case H_:
            return mtrxH__;
          case V_H:
            return mtrxV_H_;
          case VH_:
            return mtrxVH__;
        }
				return 0;
    }

	public:
    TAD_DP_TableAffine(unsigned long rows, unsigned long cols, R init)
      : TAD_DP_Table<R>(rows,cols,init) {
      checkSpaceConsumption();
      mtrxS_ = new R[this->mtrxSize_];
      mtrxV_ = new R[this->mtrxSize_];
      mtrxH_ = new R[this->mtrxSize_];
      mtrxV__ = new R[this->mtrxSize_];
      mtrxH__ = new R[this->mtrxSize_];
      mtrxV_H_ = new R[this->mtrxSize_];
      mtrxVH__ = new R[this->mtrxSize_];
			// TODO quicker way to init array?
			//std::cout << "filling with " << init << std::endl;
			std::fill( mtrxS_, mtrxS_ + this->mtrxSize_, init );
			std::fill( mtrxV_, mtrxV_ + this->mtrxSize_, init );
			std::fill( mtrxH_, mtrxH_ + this->mtrxSize_, init );
			std::fill( mtrxV__, mtrxV__ + this->mtrxSize_, init );
			std::fill( mtrxH__, mtrxH__ + this->mtrxSize_, init );
			std::fill( mtrxV_H_, mtrxV_H_ + this->mtrxSize_, init );
			std::fill( mtrxVH__, mtrxVH__ + this->mtrxSize_, init );
    }

		void checkSpaceConsumption() {
			// check for an overflow
			if (this->rows_ > ULONG_MAX / this->cols_) {
				std::cerr << "Error: Overflow in calculation matrix multiplication. Calculation terminated." << std::endl;
				exit(EXIT_FAILURE);
			}
			// maximum array size for 7 arrays is 2GB
			if ((7*this->mtrxSize_) > 2000000000 || (7*this->rows_) > 2000000000) {
        std::cerr << "Error: Maximum size of 2GB for the calculation tables exceeded due to large input data. Calculation terminated." << std::endl;
        exit(EXIT_FAILURE);
			}
		}

		inline R getMtrxVal(int table, const unsigned long i, const unsigned long j) const {
        assert(this->rowStart_[i] + j < this->mtrxSize_);
				R *mtrx = getMtrx(table);
				return mtrx[this->rowStart_[i] + j];
    }

		// TODO alg noch nicht am start
    inline void setMtrxVal(int table, const unsigned long i, const unsigned long j, const R val) {
      	assert(this->rowStart_[i] + j < this->mtrxSize_);
				R *mtrx = getMtrx(table);
				mtrx[this->rowStart_[i] + j] = val;
    }

    void print(std::ostream &s) const {
			for (int table = S; table <= VH_; table++) {
				s << table_name[table] << std::endl;
				R *mtrx = getMtrx(table);
				for (unsigned int i = 0; i < this->rows_; i++) {
					for (unsigned int j = 0; j < this->cols_; j++) {
						 s << mtrx[this->rowStart_[i] + j] << " ";
					}
					s << std::endl;
				}
				s << std::endl;
			}
			s << std::endl;
		}
};


#endif // _ALIGNMENT_TABLES_H_
