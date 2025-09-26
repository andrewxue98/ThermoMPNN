import os
from omegaconf import OmegaConf

from thermompnn.train_thermompnn import TransferModelPL


def load_model(vanilla_model_path, thermompnn_model_path):
    cfg = OmegaConf.create(
        {
            "training": {
                "num_workers": 8,
                "learn_rate": 0.001,
                "epochs": 100,
                "lr_schedule": True,
            },
            "model": {
                "hidden_dims": [64, 32],
                "subtract_mut": True,
                "num_final_layers": 2,
                "freeze_weights": True,
                "load_pretrained": True,
                "lightattn": True,
                "lr_schedule": True,
            },
            "platform": {
                "accel": "gpu",
                "thermompnn_dir": "",
                "vanilla_model_path": vanilla_model_path,
                "thermompnn_model_path": thermompnn_model_path,
            },
        }
    )

    return TransferModelPL.load_from_checkpoint(
        cfg.platform.thermompnn_model_path, cfg=cfg
    ).model
