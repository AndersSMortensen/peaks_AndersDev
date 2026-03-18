from peaks.core.fileIO.base_data_classes.base_data_class import BaseDataLoader
from peaks.core.metadata.base_metadata_models import OpticsMetadataModel

class BaseOpticsDataLoader(BaseDataLoader):
    """Baseclassfordataloders for systems with nano-focusing optics.

    Subclasses should define the `_load_optics_metadata` method to return a dictionary of relevant metadata
    values with keys of the form `temperature_item` where `item` is the names in the `_optics_attributes` list,
    i.e. is given in :class:`peaks` convention. This method should return values as :class:`pint.Quantity` objects
    where possible to ensure units are appropriately captured and propagated. Alternatively, the main `_load_metadata`
    method can be overwritten to return the full metadata dictionary, including manipulator metadata.

    Subclasses should add any additional temperature attributes via the `_add_temperature_attributes` class variable,
     providing a list of additional attributes.

    See Also
    --------
    BaseDataLoader
    BaseDataLoader._load_metadata
    """

    #Define class variables
    _loc_name = "Default Optics"
    _optics_attributes = ["opt_x1","opt_x2","opt_x3"]
    _metadata_parsers = [
        "_parse_optics_metadata"
    ]

    #Properties to access class variables
    @property
    def tempetrature_attributes(self):
        "Return the optics attributes."
        return self._optics_attributes
    
    @classmethod
    def _parse_optics_metadata(cls, metadata_dict):
        """Parse metadataspecific to the optics data."""

        #Build and populate the optics metadata model.
        optics_metadata = OpticsMetadataModel(
            opt_x1=metadata_dict.get("optics_opt_x1"),
            opt_x2=metadata_dict.get("optics_opt_x1"),
            opt_x3=metadata_dict.get("optics_opt_x3"),
        )

        metadata_to_warn_if_missing = (
            f"optics_{attribute}"
            for attribute in cls._optics_attributes
        )

        return {"_optics":optics_metadata}, metadata_to_warn_if_missing