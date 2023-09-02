module UnitfulExt

import COESA:
    COESA,
    atmosphere,
    altitude,
    mean_molecular_weight,
    temperature,
    pressure,
    density,
    speed_of_sound,
    dynamic_viscosity

isdefined(Base, :get_extension) ? (using Unitful) : (using ..Unitful)

struct UnitfulState
    state::COESA.State
end

altitude(s::UnitfulState) = altitude(s.state) * u"m"
mean_molecular_weight(s::UnitfulState) = mean_molecular_weight(s.state) * u"kg/kmol"
temperature(s::UnitfulState) = temperature(s.state) * u"K"
pressure(s::UnitfulState) = pressure(s.state) * u"Pa"
density(s::UnitfulState) = density(s.state) * u"kg/m^3"
speed_of_sound(s::UnitfulState) = speed_of_sound(s.state) * u"m/s"
dynamic_viscosity(s::UnitfulState) = dynamic_viscosity(s.state) * u"N*s/m^2"

function atmosphere(Z::Quantity)
    Z′ = ustrip(u"m", Z)
    state = atmosphere(Z′)
    UnitfulState(state)
end

end
