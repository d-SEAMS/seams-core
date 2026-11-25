{ fetchurl, system }:

let
  wheel = url: sha256: fetchurl { inherit url sha256; };
  macQt = {
    pyside = wheel
      "https://files.pythonhosted.org/packages/e1/90/d528f8dd3f9f20712bc719e10c698d27c998e1dd76ed3f33b6fb606d0972/pyside6_essentials-6.12.0-cp310-abi3-macosx_14_0_universal2.whl"
      "0867b709a724db28b161227c7d6ecb5ed10cc09226a53542417ac6dee1308239";
    shiboken = wheel
      "https://files.pythonhosted.org/packages/b5/45/11536fba5f50a57671eb14441c6098096f3e57074e54c9d7d66039a5831e/shiboken6-6.12.0-cp310-abi3-macosx_14_0_universal2.whl"
      "47105e05baf57d35453e07d240bbe59cb4b8d9565377ea5134ff9fb37d4f5c02";
  };
  wheels = {
    x86_64-linux = {
      ovito = wheel
        "https://files.pythonhosted.org/packages/67/99/09355940b7be501f016e481e69c67aa795aac706dd749b0b34ca6f6be482/ovito-3.17.2-cp314-cp314-manylinux_2_34_x86_64.whl"
        "1j6zplnpzyn7s87h79cgqkilfpp4ganbb5xv86vwcbz4sdwl8sr3";
      pyside = wheel
        "https://files.pythonhosted.org/packages/7b/3c/7597d8dcb36aa8315427635ec6ec9d850ffbad4ce78b9118ac854889f2d3/pyside6_essentials-6.12.0-cp310-abi3-manylinux_2_34_x86_64.whl"
        "1n9v7n4nkkjllidcrya0h0ymzzlyv91q4ia8w1d26jdg3pjbbx0n";
      shiboken = wheel
        "https://files.pythonhosted.org/packages/49/c2/a04fa003b4360e99d8e3fc2829abacf37396ace63ea1cd616402209ec7c6/shiboken6-6.12.0-cp310-abi3-manylinux_2_34_x86_64.whl"
        "0h329n41j6fnzqmmlfg9871aprgdvnnhhf7nqysl8049n1dmdpiv";
    };
    aarch64-linux = {
      ovito = wheel
        "https://files.pythonhosted.org/packages/38/06/41be9c68650da0c4732f9adc45ee5693679bd4a4a523c826d93a1fedc022/ovito-3.17.2-cp314-cp314-manylinux_2_39_aarch64.whl"
        "dcc3d836df058bba7b648447557dc31dff3e0fb76230102b7ca10d7b750df7da";
      pyside = wheel
        "https://files.pythonhosted.org/packages/5e/4e/f4fd4c016719be832ec49d3c87940459e87fd666179a0a199a3998d07916/pyside6_essentials-6.12.0-cp310-abi3-manylinux_2_39_aarch64.whl"
        "08c5841063fc1df69b7fa1ee405f539e66576a83fc48004be3e2c5e35ccbfa69";
      shiboken = wheel
        "https://files.pythonhosted.org/packages/97/ca/a4415e6218e0dff764fe2b2d30db2fc8f499a08ace07634232ade3f73adc/shiboken6-6.12.0-cp310-abi3-manylinux_2_39_aarch64.whl"
        "ff72a72b3277b418902d562b987b33a6f6818bed9fe40a49d634d09bc67f620e";
    };
    aarch64-darwin = macQt // {
      ovito = wheel
        "https://files.pythonhosted.org/packages/0b/39/3e936b6ae737b375b9b3b73ed6636c58f42a449d968c2e5ae26770649d97/ovito-3.17.2-cp314-cp314-macosx_14_0_arm64.whl"
        "95002125c870f2c593eb075a77528e624770a6ea28910e40811157d835ac38aa";
    };
    x86_64-darwin = macQt // {
      # OVITO 3.12 is the last release with an Intel macOS wheel.
      ovito = wheel
        "https://files.pythonhosted.org/packages/a5/62/8f4b630e2b14cfc6e2395f5c8ef1c2b829ffda756025302d28cdcb33a25d/ovito-3.12.0-cp313-cp313-macosx_11_0_x86_64.whl"
        "cd321fe555d7e0f2b2f91b6fec97c25f8597f158d12faf62f47f443cbac68563";
    };
  };
in
wheels.${system}
