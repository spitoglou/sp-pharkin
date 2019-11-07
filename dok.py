import sp_pharkin.clearance as sppkcl
import sp_pharkin as sppk


if __name__ == '__main__':

    v = sppkcl.average_clearance_weight(
        average_clearance='0.04L/(hour*kilogram)',
        weight='75 kilogram'
    )
    print(v)

    print('4.1')
    v = sppk.half_life_k(
        K='0.0131 /hour',
    )
    print(v)

    print('4.2')
    v = sppkcl.clearance_flow_extraction_rate(
        Q='1.2L/min',
        E=0.04,
        output_unit='L/hour'
    )
    print(v)

    print('4.3')
    v = sppk.half_life_k(
        half_life='12hour',
        decimals=3
    )
    print(v)

    print('4.4')
    v = sppkcl.clearance_elimination_rate_constant_volume(
        volume='25L',
        K='0.03/hour',
        output_unit='mL/hour'
    )
    print(v)

    print('4.5')
    k = sppkcl.clearance_elimination_rate_constant_volume(
        volume='175L',
        clearance='0.43L/min',
        decimals=10,
        output_unit='1/hour'
    )
    print(k)
    v = sppk.half_life_k(
        K=k[3],
        output_unit='hour'
    )
    print(v)

    print('4.6')
    k = sppkcl.clearance_elimination_rate_constant_volume(
        volume='12L',
        clearance='12.5mL/min',
        decimals=4,
        output_unit='1/hour'
    )
    print(k)

    print('4.7')
    k = sppkcl.average_clearance_weight(
        average_clearance='1.5mL/min/kilogram',
        weight='85kilogram',
        output_unit='L/hour'
    )
    print(k)
