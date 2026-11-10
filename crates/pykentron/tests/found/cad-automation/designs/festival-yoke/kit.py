"""Builders for the modules in the active kit (see layout.py and config.toml [kits])."""
import layout
import module
import power

BUILDERS = {
    "cup": module.cup_module, "phone": module.phone_module, "dock": module.dock_module, "fan": module.fan_module,
    "pod": power.pod_module, "id": module.id_module, "tube": module.tube_module,
    "blank": module.blank_module,
}


def modules():
    """name -> shape for every module in the active kit, at its tuned angle."""
    return {name: BUILDERS[name]() for name in layout.present()}
