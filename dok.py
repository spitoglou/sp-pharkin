import sp_pharkin.clearance as sppkcl
import sp_pharkin.expo as expo
# import sp_pharkin as sppk


if __name__ == '__main__':

    v = sppkcl.average_clearance_weight(
        average_clearance='0.04L/(hour*kilogram)',
        weight='75 kilogram'
    )
    print(v)

    expo.c_t('k', c_0=100, c_t=50, t=20)
