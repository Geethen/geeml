# import rasterio as rio
# import math
# import logging
# from concurrent.futures import ThreadPoolExecutor, as_completed
# import concurrent.futures
# import threading
# import torch
# from tqdm.auto import tqdm

# def denormalise_output(image):
#     image = torch.nan_to_num(image, nan=0.0, posinf= 0.0, neginf= 0.0)
#     min_i = torch.full((image.shape[0], image.shape[1]), target_min[0], device = device)
#     max_i = torch.full((image.shape[0], image.shape[1]), target_max[0], device = device)
#     denorm_img = (image * (max_i - min_i)) + min_i
#     return denorm_img

# def inference(infile, model, outfile, patchSize, num_workers=4):
#     """
#     Run inference using model on infile block-by-block and write to a new file (outfile). 
#     In the case, that the infile image width/height is not exactly divisible by 32, padding
#     is added for inference and removed prior to the outfile being saved.
    
#     Args:
#         infile (string): Path to input image/covariates
#         model (pth file): Loaded trained model/checkpoint
#         outfile (string): Path to save predicted image
#         patchSize (int): Must be a multiple of 32. Size independent of model input size.
#         num_workers (int): Num of workers to parralelise across
        
#     Returns:
#         A tif saved to the outfile destination
        
#     """

#     with rio.open(infile) as src:
        
#         logger = logging.getLogger(__name__)

#         # Create a destination dataset based on source params. The
#         # destination will be tiled, and we'll process the tiles
#         # concurrently.
#         profile = src.profile
#         profile.update(blockxsize= patchSize, blockysize= patchSize, tiled=True, count=1)

#         with rio.open(Path(outfile), "w", **profile) as dst:
#             windows = [window for ij, window in dst.block_windows()]

#             # use a lock to protect the DatasetReader/Writer
#             read_lock = threading.Lock()
#             write_lock = threading.Lock()

#             def process(window):
#                 with read_lock:
#                     src_array = src.read(window=window)#nbands, nrows, ncols(4, h, w)
#                     w, h = src_array.shape[1], src_array.shape[2]
#                     image = normalise_input(torch.from_numpy(src_array))
#                     image = torch.from_numpy(np.expand_dims(image, axis=0)).to(device, dtype=torch.float)#(1, h, w, 4)
#                     hpad = math.ceil(h/32)*32-h
#                     wpad = math.ceil(w/32)*32-w
#                     transform = transforms.Pad((0, 0, hpad, wpad))
#                     # add padding to image
#                     image = transform(image)
#                     model.eval()
#                     output = model(image)#(1,1,h,w)
#                     #denormalize output
#                     output = denormalise_output(output.squeeze())
#                     # remove padding
#                     result = output[0:w, 0:h].detach().cpu()
#                     # plt.imshow(result.numpy())

#                 with write_lock:
#                     dst.write(result, 1, window=window)

#             # We map the process() function over the list of
#             # windows.
#             with tqdm(total=len(windows), desc = os.path.basename(outfile)) as pbar:
#                 with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
#                     futures = {executor.submit(process, window): window for window in windows}
                    
#                     try:
#                         for future in concurrent.futures.as_completed(futures):
#                             future.result()
#                             pbar.update(1)
                                    
#                     except Exception as ex:
#                         logger.info('Cancelling...')
#                         executor.shutdown(wait=False, cancel_futures=True)
#                         raise ex        